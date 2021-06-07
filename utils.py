from typing import *

import numpy as np

import err_corr_code
from constants import *




def to_bytestring(data: str):
    if not isinstance(data, bytes):
        data = data.encode("utf-8")
    return data

class RSBlock(object):

    def __init__(self, total_count, data_count):
        self.total_count = total_count
        self.data_count = data_count


def rs_blocks(version, error_correction):
    if error_correction not in ERROR_CORR_LEVEL:  # pragma: no cover
        raise Exception(
            "bad rs block @ version: %s / error_correction: %s" %
            (version, error_correction))
    offset = ERROR_CORR_LEVEL[error_correction]
    rs_block = RS_BLOCK_TABLE[(version - 1) * 4 + offset]

    blocks = []

    for i in range(0, len(rs_block), 3):
        count, total_count, data_count = rs_block[i:i + 3]
        for j in range(count):
            blocks.append(RSBlock(total_count, data_count))

    return blocks

class BitBuffer:
    """a library to store one value (0-255) bit by bit
    """
    def __init__(self) -> None:
        self.buffer = []
        self.length = 0
        
    def __len__(self):
        return self.length
    
    def get(self, index):
        
        import math
        buf_index = math.floor(index / 8)
        return ((self.buffer[buf_index] >> (7 - index % 8)) & 1) == 1
    def put(self, num, length):
        for i in range(length):
            self.put_bit(((num >> (length - i - 1)) & 1) == 1)
    
    def put_bit(self, bit):
        buf_index = self.length // 8
        if len(self.buffer) <= buf_index:
            self.buffer.append(0)
        if bit:
            self.buffer[buf_index] |= (0x80 >> (self.length % 8))
        self.length += 1

_helper_func = lambda x: x.data_count
VERSION_ERRCORR_BITS = [[0] + [8 * sum(map(_helper_func, rs_blocks(version, corr_level))) for version in range(1, 41)] for corr_level in ["M", "L", "H", "Q"]]


class QrData:
    def __init__(self, data, check_data=True) -> None:
        if check_data:
            data = to_bytestring(data)
        self.data = data
        
    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return repr(self.data)
    
    def write(self, buffer: BitBuffer):
        for c in self.data:
            buffer.put(c, 8)




def BCH_digit(data):
    """count bits in binary representation of data

    Args:
        data ([int]): [description]

    Returns:
        [int]: [digits]
    """
    digit = 0
    while data != 0:
        digit += 1
        data >>= 1
    return digit

     
def BCH_format_code(data):
    d = data << 10
    while BCH_digit(d) >= BCH_digit(G15):  # polynomials divide operation
        d ^= (G15 << (BCH_digit(d) - BCH_digit(G15)))  # binary modulus: minus must be treated as plus
    return ((data << 10) | d) ^ G15_MASK  # follosing annex C in QR code Spec.
    
def BCH_version_code(data):
    d = data << 12
    while BCH_digit(d) >= BCH_digit(G18):
        d ^= (G18 << (BCH_digit(d) - BCH_digit(G18)))
    return (data << 12) | d

def mask_function(pattern):
    assert 0 <= pattern <= 7, "Bad mask pattern: {}, expected 0-7".format(pattern)
    if pattern == 0:
        return lambda i, j: (i + j) % 2 == 0
    elif pattern == 1:
        return lambda i, j: i % 2 == 0
    elif pattern == 2:
        return lambda i, j: j % 3 == 0
    elif pattern == 3:
        return lambda i, j: (i + j) % 3 == 0
    elif pattern == 4:
        return lambda i, j: ((i // 2) + (j // 3)) % 2 == 0
    elif pattern == 5:
        return lambda i, j: (i * j) % 2 + (i * j) % 3 == 0
    elif pattern == 6:
        return lambda i, j: ((i * j) % 2 + (i * j) % 3) % 2 == 0
    elif pattern == 7:
        return lambda i, j: ((i * j) % 3 + (i + j) % 2) % 2 == 0

def count_lost_norm1(matrix):
    """calculate 1st lost points: Adjacent modules in row/column in same color: 3 + (length - 5)

    Args:
        matrix ([type]): [description]

    Returns:
        [int]: [description]
    """
    n, _ = matrix.shape
    
    lost_point = 0
    for row in range(n):
        previous_color = matrix[row, 0]
        length = 0
        for col in range(n):
            if matrix[row, col] != previous_color:
                if length >= 5:
                    lost_point += length - 2
                length = 1
                previous_color = matrix[row, col]
            else:
                length += 1
        if length >= 5:
            lost_point += length - 2
        
    for col in range(n):
        previous_color = matrix[0, col]
        length = 0
        for row in range(n):
            if matrix[row, col] != previous_color:
                if length >= 5:
                    lost_point += length - 2
                length = 1
                previous_color = matrix[row, col]
            else:
                length += 1
        if length >= 5:
            lost_point += length - 2
    return lost_point
    
        
def count_lost_norm2(matrix):
    """calculate 2nd lost points: Block of modules in same color: 3 * (m - 1) * (n - 1)

    Args:
        matrix ([type]): [description]

    Returns:
        [type]: [description]
    """
    n, _ = matrix.shape
    lost_point = 0
    for row in range(n - 1):
        this_row = matrix[row]
        next_row = matrix[row + 1]
        # use iterator to reduce time
        col_iterator = iter(range(n - 1))
        for col in col_iterator:
            top_right = this_row[col + 1]
            if top_right != next_row[col + 1]:
                next(col_iterator, None)  # iterator next, then after the loop, iterator next again
            elif top_right != this_row[col]:
                continue
            elif top_right != next_row[col]:
                continue
            else:
                lost_point += 3  # N2 = 3, each col will increase three lost points
    return lost_point
                

def count_lost_norm3(matrix):
    """calculate 3rd lost points: 1:1:3:1:1 ratio (dark:light:dark:light:dark) pattern
    in row/column. N3 = 40

    Args:
        matrix ([type]): [description]

    Returns:
        [type]: [description]
    """
    # only two types: 10111010000 and 00001011101
    # reference from ISOIEC
    n, _ = matrix.shape
    lost_point = 0
    normal_range = range(n)
    short10_range = range(n - 10)
    for row in normal_range:
        short10_iter = iter(short10_range)
        for col in short10_iter:
            if (
                # 10111010000
                matrix[row, col] and not matrix[row, col + 1] and matrix[row, col + 2]
                and matrix[row, col + 3] and matrix[row, col + 4] and not matrix[row, col + 5]
                and matrix[row, col + 6] and not matrix[row, col + 7] and not matrix[row, col + 8]
                and not matrix[row, col + 9] and not matrix[row, col + 10]
            ) or (
                # 00001011101
                not matrix[row, col] and not matrix[row, col + 1] and not matrix[row, col + 2]
                and not matrix[row, col + 3] and matrix[row, col + 4] and not matrix[row, col + 5]
                and matrix[row, col + 6] and matrix[row, col + 7] and matrix[row, col + 8]
                and not matrix[row, col + 9] and matrix[row, col + 10]
            ):
                lost_point += 40
            
            
            if matrix[row, col + 10]:
                # horspool algorithm for prefix match
                next(short10_iter, None) # give a default value to avoid exception
    for col in normal_range:
        short10_iter = iter(short10_range)
        for row in short10_iter:
            if (
                matrix[row, col] and not matrix[row + 1, col] and matrix[row + 2, col]
                and matrix[row + 3, col] and matrix[row + 4, col] and not matrix[row + 5, col]
                and matrix[row + 6, col] and not matrix[row + 7, col] and not matrix[row + 8, col]
                and not matrix[row + 9, col] and not matrix[row + 10, col]
            ) or (
                # 00001011101
                not matrix[row, col] and not matrix[row + 1, col] and not matrix[row + 2, col]
                and not matrix[row + 3, col] and matrix[row + 4, col] and not matrix[row + 5, col]
                and matrix[row + 6, col] and matrix[row + 7, col] and matrix[row + 8, col]
                and not matrix[row + 9, col] and matrix[row + 10, col]
            ):
                lost_point += 40
            if matrix[row + 10, col]:
                next(short10_iter, None)
    return lost_point 

def count_lost_norm4(matrix):
    """calculate 4th lost points: Proportion of dark modules in entire symbol: 
    50 + (5 + k) or 50 - (5 + k), return k * 10

    Args:
        matrix ([type]): [description]

    Returns:
        [int]: [description]
    """
    dark_sum = np.sum(matrix)
    modules_num = matrix.size
    dark_ratio = dark_sum / modules_num
    k = abs((dark_ratio * 100 - 50)) / 5
    return int(k) * 10    

def count_lost(matrix):
    """count lost points in filled matrix. Modules have been masked by some mask pattern

    Args:
        matrix ([np.ndarray]): [description]

    Returns:
        [int]: [how many points are lost]
    """
    return count_lost_norm1(matrix) + count_lost_norm2(matrix) + count_lost_norm3(matrix) + count_lost_norm4(matrix)
    
    
def generate_data(version, err_corr, data_list):
    """first add terminator, then rearange by 8 bits with padding. Finally, add error correction codes
    and compose to an entity

    Args:
        version ([int]): [description]
        err_corr ([int]): [1--L, 0--M, 3--Q, 2--H]
        data_list ([list[bytes]]): [description]

    Returns:
        [type]: [description]
    """
    each_binary_term_len = BITS_CHAR_COUNT_INDICATOR(version)
    buffer = BitBuffer()
    for data in data_list:
        buffer.put(MODE_INDICATOR, 4)
        buffer.put(len(data), each_binary_term_len)
        
        data.write(buffer)
    
    filled_err_corr = "L" if err_corr == 1 else "M" if err_corr == 0 else "Q" if err_corr == 3 else "L"
    blocks = rs_blocks(version, filled_err_corr)
    bit_limit = sum([8 * block.data_count for block in blocks])
    
    assert bit_limit >= len(buffer), "Code length overflow. Data size {} > size available ({})".format(len(buffer), bit_limit)
    
    # fill terminator
    for i in range(min(bit_limit - len(buffer), 4)):
        buffer.put_bit(False)
    # rearange by 8 bit
    remaining = len(buffer) % 8
    if remaining > 0:
        for i in range(8 - remaining):
            buffer.put_bit(False)
    # padding
    padding_bytes = (bit_limit - len(buffer)) // 8

    for i in range(padding_bytes):
        buffer.put(PAD1 if i % 2 == 0 else PAD2, 8)
    
    return generate_bytes(buffer, blocks)

def generate_bytes(buffer: BitBuffer, rs_blocks: List[RSBlock]):
    """part of this function is referred to https://github.com/lincolnloop/python-qrcode
    a normal way to generative error correction codeword using Polynomials

    Args:
        buffer (BitBuffer): [description]
        rs_blocks (List[RSBlock]): [description]

    Returns:
        [List[int]]: [data codeword in 1-D array]
    """
    offset = 0  # offset to divide blocks for data codewords
    max_data_num = 0
    max_err_num = 0
    
    block_num = len(rs_blocks)
    encoded_data = [0] * block_num
    encoded_err = [0] * block_num
    for block in range(block_num):
        block_data_num = rs_blocks[block].data_count
        block_err_num = rs_blocks[block].total_count - block_data_num
        
        max_data_num = max(max_data_num, block_data_num)
        max_err_num = max(max_err_num, block_err_num)
        
        encoded_data[block] = [0] * block_data_num
        
        for i in range(len(encoded_data[block])):
            encoded_data[block][i] = 255 & buffer.buffer[i + offset]  # avoid the power of g(x) exceeds 255
        
        offset += block_data_num
        
        # following part exceeds my ability. I understande the generation method and process
        # but cannot complete the code by myself
        # reference from QRCoda Spec.
        if block_err_num in gen_poly_coff:
            gen_poly = err_corr_code.Polynomial(gen_poly_coff[block_err_num], 0)
        else:
            gen_poly = err_corr_code.Polynomial([1], 0)
            for i in range(block_err_num):
                gen_poly = gen_poly * err_corr_code.Polynomial([1, err_corr_code.gexp(i)], 0)
        # after getting generative polynomials for certain number of error correction codewords
        # long divide method to get modulus
        message_poly = err_corr_code.Polynomial(encoded_data[block], len(gen_poly) - 1)
        
        raw_err_corr_code = message_poly % gen_poly
        encoded_err[block] = [0] * (len(gen_poly) - 1)
        cur_len = len(encoded_err[block])
        for i in range(cur_len):
            raw_idx = i + len(raw_err_corr_code) - cur_len
            if raw_idx >= 0:
                encoded_err[block][i] = raw_err_corr_code[raw_idx]
            else:
                encoded_err[block][i] = 0
        # till now, error correction codeword is generated
    # composite to the whole data codeword
    
    total_code_num = sum([block.total_count for block in rs_blocks])
    data = [0] * total_code_num
    total_idx = 0
    # fill data codeword
    for i in range(max_data_num):
        for j in range(len(rs_blocks)):
            # for each block, place first data
            # then for next data
            if i < len(encoded_data[j]):
                data[total_idx] = encoded_data[j][i]
                total_idx += 1
    
    # fill err corr codeword
    for i in range(max_err_num):
        for j in range(len(rs_blocks)):
            if i < len(encoded_err[j]):
                data[total_idx] = encoded_err[j][i]
                total_idx += 1
    return data



if __name__ == "__main__":
    print(rs_blocks(5, "H"))
