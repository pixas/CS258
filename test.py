import math

class BitBuffer(object):

    def __init__(self):
        self.buffer = []
        self.length = 0 

    def __repr__(self):
        return ".".join([str(n) for n in self.buffer])

    def get(self, index):
        buf_index = math.floor(index / 8)
        return ((self.buffer[buf_index] >> (7 - index % 8)) & 1) == 1

    def put(self, num, length):
        for i in range(length):
            self.put_bit(((num >> (length - i - 1)) & 1) == 1)

    def __len__(self):
        return self.length

    def put_bit(self, bit):
        buf_index = self.length // 8
        if len(self.buffer) <= buf_index:
            self.buffer.append(0)
        if bit:
            self.buffer[buf_index] |= (0x80 >> (self.length % 8))
        self.length += 1
        
        
        

if __name__ == "__main__":
    import qrcode
    q = qrcode.QRCode()
    q.add_data("I have a dream")
    print(q.get_matrix())