EXP_TABLE = list(range(256))

LOG_TABLE = list(range(256))

for i in range(8):
    EXP_TABLE[i] = 1 << i

for i in range(8, 256):
    EXP_TABLE[i] = (
        EXP_TABLE[i - 4] ^ EXP_TABLE[i - 5] ^ EXP_TABLE[i - 6] ^
        EXP_TABLE[i - 8])

for i in range(255):
    LOG_TABLE[EXP_TABLE[i]] = i
    

def glog(n):
    """from an value infer the power of item in polynomial

    Args:
        n ([int]): [coefficient]

    Raises:
        ValueError: [coefficient must larget or equal to 1]

    Returns:
        [int]: [range from 0 to 255]
    """
    if n < 1:  # pragma: no cover
        raise ValueError("glog(%s)" % n)
    return LOG_TABLE[n]


def gexp(n):
    """get R field number for nth power coefficient

    Args:
        n ([int]): [alpha's power. may larger than 255]

    Returns:
        [type]: [description]
    """
    return EXP_TABLE[n % 255]


class Polynomial(object):

    def __init__(self, num, shift):
        if not num:  # pragma: no cover
            raise Exception("%s/%s" % (len(num), shift))

        for offset in range(len(num)):
            if num[offset] != 0:
                break
        else:
            offset += 1
        # self.num means coefficients
        # shift means error correction code's number
        self.num = num[offset:] + [0] * shift

    def __getitem__(self, index):
        return self.num[index]

    def __iter__(self):
        return iter(self.num)

    def __len__(self):
        return len(self.num)

    def __mul__(self, other):
        num = [0] * (len(self) + len(other) - 1)

        for i, item in enumerate(self):
            for j, other_item in enumerate(other):
                num[i + j] ^= gexp(glog(item) + glog(other_item))
                # multiply in GF means add(xor) corresponding coefficients

        return Polynomial(num, 0)

    def __mod__(self, other):
        difference = len(self) - len(other)
        if difference < 0:
            return self
        # highest order item coefficient
        ratio = glog(self[0]) - glog(other[0])

        num = [
            item ^ gexp(glog(other_item) + ratio)  # multiply in R means add in GF
            for item, other_item in zip(self, other)]
        if difference:  # if self's order > other's order in highest item
            num.extend(self[-difference:])

        # recursive call
        return Polynomial(num, 0) % other