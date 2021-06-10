import os
from bisect import bisect_left

import numpy as np
from matplotlib import pyplot as plt

from constants import *
from utils import *
import cv2

class QrCode:
    def __init__(self, version=None, err_corr="M",
                 border=4, mask_pattern=None) -> None:

        self.version = version and int(version)
        self.err_corr = 1 if err_corr == "L" else 0 if err_corr == "M" else 3 if err_corr == "Q" else 2
        self.mask_pattern = mask_pattern
        # self.excel = pd.read_excel("GenerationConstants.xlsx", "Number of Symbols")
        # with open("encoding table.pkl", 'rb') as file:
        #     self.encode_table = pickle.load(file)
        self.border = border
        self.reset()

    def reset(self):
        self.modules = None
        self.modules_num = 0
        self.data_cache = None
        self.data_list = []

    def add_data(self, data):
        # TODO yet not apply to all circumstances
        # will complete other circumstances after the final week
        if isinstance(data, QrData):
            self.data_list.append(data)
        else:
            # pre disposal
            # now self.data_list stores one information
            # if applied to all kinds of data, this list contains more than one elements for 
            # optimal splits by the method in annex of QR Spec.
            self.data_list.append(QrData(data))
        self.data_cache = None

    def make(self, fit=True):
        if fit or (self.version is None):
            self.best_fit(self.version)
        if self.mask_pattern is None:
            self.generate(False, self.best_mask_pattern())
        else:
            self.generate(False, self.mask_pattern)

    def best_fit(self, start):
        if start is None:
            start = 1

        # convert to binary sequences

        each_binary_term_len = BITS_CHAR_COUNT_INDICATOR(start)
        buffer = BitBuffer()
        for data in self.data_list:
            buffer.put(MODE_INDICATOR, 4)
            buffer.put(len(data), each_binary_term_len)
            data.write(buffer)

        now_bits = len(buffer)

        self.version = bisect_left(
            VERSION_ERRCORR_BITS[self.err_corr], now_bits, start)
        assert 1 <= self.version <= 40, "Overflow error! expects version in [1, 40] range, got {}".format(
            self.version)
        # test this version match or not
        if each_binary_term_len != BITS_CHAR_COUNT_INDICATOR(self.version):
            self.best_fit(self.version)
        return self.version

    def best_mask_pattern(self):
        mask_pattern = 0
        min_lost = 0

        for i in range(8):
            self.generate(True, mask_pattern)

            cur_mask_lost = count_lost(self.modules)
            if i == 0 or cur_mask_lost < min_lost:
                min_lost = cur_mask_lost
                mask_pattern = i
        return mask_pattern

    def fill_position_detection_pattern(self, row, col):
        for r in range(-1, 8):
            if r + row < 0 or r + row >= self.modules_num:
                continue
            for c in range(-1, 8):
                if c + col < 0 or c + col >= self.modules_num:
                    continue
                if (0 <= r <= 6 and (c == 0 or c == 6)) \
                        or (0 <= c <= 6 and (r == 0 or r == 6)) \
                        or (2 <= r <= 4 and 2 <= c <= 4):
                    self.modules[row + r, col + c] = True
                else:
                    self.modules[row + r, col + c] = False

    def fill_alignment_pattern(self):
        positions = PATTERN_POSITION_TABLE[self.version - 1]
        length = len(positions)
        for i in range(length):
            for j in range(length):
                row = positions[i]
                col = positions[j]

                if self.modules[row, col] is not None:
                    continue
                for r in range(-2, 3):
                    for c in range(-2, 3):
                        if (r == -2 or r == 2 or c == -2 or c == 2 or (c == 0 and r == 0)):
                            self.modules[row + r, col + c] = True
                        else:
                            self.modules[row + r, col + c] = False
        return

    def fill_timing_pattern(self):
        for col in range(8, self.modules_num - 8):
            if self.modules[6, col] is not None:
                continue
            self.modules[6, col] = (col % 2 == 0)
        for row in range(8, self.modules_num - 8):
            if self.modules[row, 6] is not None:
                continue
            self.modules[row, 6] = (row % 2 == 0)

        return

    def fill_format_information(self, test, mask_pattern):
        format_info = (self.err_corr << 3) | mask_pattern
        data_plus_BCH = BCH_format_code(format_info)

        for row in range(0, 15):
            filled_val = (not test and ((data_plus_BCH >> row) & 1) == 1)
            if row < 6:
                self.modules[row, 8] = filled_val
            elif row < 8:
                self.modules[row + 1, 8] = filled_val
            else:
                self.modules[self.modules_num - 15 + row, 8] = filled_val

        for col in range(0, 15):
            filled_val = (not test and ((data_plus_BCH >> col) & 1) == 1)
            if col < 8:
                self.modules[8, self.modules_num - col - 1] = filled_val
            elif col < 9:
                self.modules[8, 15 - col] = filled_val
            else:
                self.modules[8, 15 - col - 1] = filled_val
        # forever dark
        self.modules[self.modules_num - 8, 8] = not test
        return

    def fill_version_information(self, test):
        """filled the qr code. This function is only called when self.version >= 7

        Args:
            test ([bool]): [description]
        """
        data_plus_BCH = BCH_version_code(self.version)
        for i in range(18):
            filled_val = (not test and ((data_plus_BCH >> i) & 1) == 1)
            self.modules[i // 3, self.modules_num - 11 + i % 3] = filled_val
        for i in range(18):
            filled_val = (not test and ((data_plus_BCH >> i) & 1) == 1)
            self.modules[self.modules_num - 11 + i % 3, i // 3] = filled_val
        return

    def fill_data_err_corr_code(self, data, mask_pattern):
        """fill the modules by Z shape. Remember: although some modules are filled, sequence can stll cover
        these grids but do nothing. This way can reduce code complexity

        Args:
            data ([type]): [description]
            mask_pattern ([type]): [description]
        """
        mask_func = mask_function(mask_pattern)
        data_size = len(data)
        increment = -1
        row = self.modules_num - 1
        bitIndex = 7
        byteIndex = 0
        for col in range(self.modules_num - 1, 0, -2):
            if col <= 6:
                col -= 1
            col_range = [col, col - 1]
            while True:
                for c in col_range:
                    if self.modules[row][c] is None:
                        # if filled by other information, skip
                        filled_val = False
                        if byteIndex < data_size:
                            filled_val = (
                                (data[byteIndex] >> bitIndex) & 1) == 1
                        if mask_func(row, c):
                            filled_val = not filled_val
                        self.modules[row, c] = filled_val
                        bitIndex -= 1

                        if bitIndex < 0:
                            byteIndex += 1
                            bitIndex = 7
                row += increment

                if row < 0 or row >= self.modules_num:
                    row -= increment
                    increment = -increment
                    break

        return

    def generate(self, test, mask_pattern):
        """generate qr image for mask or final image.

        Args:
            test ([bool]): [test True for mask, False for original qr code]
            mask_pattern ([int]): [0 - 7 designates the method to mask]
        """
        self.modules_num = self.version * 4 + 17
        self.modules = np.array(
            [[None for i in range(self.modules_num)] for j in range(self.modules_num)])

        # fill three position detection patterns
        self.fill_position_detection_pattern(0, 0)
        self.fill_position_detection_pattern(0, self.modules_num - 7)
        self.fill_position_detection_pattern(self.modules_num - 7, 0)
        self.fill_alignment_pattern()
        self.fill_timing_pattern()
        self.fill_format_information(test, mask_pattern)

        if self.version >= 7:
            self.fill_version_information(test)

        if self.data_cache is None:
            self.data_cache = generate_data(
                self.version, self.err_corr, self.data_list)
        self.fill_data_err_corr_code(self.data_cache, mask_pattern)

    def display(self, save_dir=None, name=None):

        array = np.pad(self.modules, ((self.border, self.border), (self.border, self.border)),
                       "constant", constant_values=((False, False), (False, False)))

        array = np.array(~array, int)

        # print(array)
        

        fig, ax = plt.subplots()
        ax.imshow(array, "gray")
        ax.set_xticks([])
        ax.set_yticks([])
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        if save_dir is not None:
            plt.savefig(os.path.join(save_dir, name), dpi=600)
        plt.show()
        
    def insert_logo(self, logo):

        pixels = np.tile(np.expand_dims(np.where(self.modules, 0, 255), -1), [1, 1, 3])
        size = self.modules.shape[0]
        if isinstance(logo, np.ndarray):
            emb_img = logo
        elif isinstance(logo, str):
            emb_img = cv2.imread(os.path.join(os.path.dirname(os.path.abspath(__file__)), logo), cv2.IMREAD_GRAYSCALE)
            print((emb_img == 1).any())
            
            emb_img = cv2.cvtColor(emb_img, cv2.COLOR_BGR2RGB)
        else:
            raise TypeError("logo argument expects str or np.ndarray dtype, got {}".format(type(logo)))
        
        plt.imshow(emb_img)
        plt.show()
        
        emb_width = int(min(emb_img.shape[1], size / 7 * 2))
        emb_height = int(min(emb_img.shape[0], size / 7 * 2))
        # print(emb_width, emb_height)
        emb_fill_height = (size - emb_height) // 2
        emb_fill_width = (size - emb_width) // 2
        
        emb_img = cv2.resize(emb_img, (emb_width, emb_height), interpolation=cv2.INTER_AREA)
        pixels[emb_fill_height: emb_fill_height + emb_height, emb_fill_width: emb_fill_width + emb_width, :] = emb_img
        plt.imshow(pixels)
        plt.show()
        pass

if __name__ == "__main__":
    # q = QrCode()
    import time
    q = QrCode()
    q.add_data("""My number is 519030910043""")
    begin = time.time()
    q.make()
    end = time.time()
    print("Execution time: {}".format(end - begin))
    print(q.version)
    # q.display(save_dir=os.path.dirname(os.path.abspath(__file__)), name="zen.jpg")
    q.display()
    # q.insert_logo("small.png")
    # import qrcode
    # q = qrcode.QRCode()
    # q.add_data("""Through a flowery sea of dreams they go.""")
    # matrix = q.get_matrix()
    # print(q.version)
    
    # plt.imshow(~np.array(matrix), "gray")
    # plt.show()
