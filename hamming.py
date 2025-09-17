## hamming distance is number of characters that disagree

## http://code.activestate.com/recipes/499304-hamming-distance/
import operator
def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

## here's a more straightforward implementation, not as fast
def hamming1(str1,str2):
    assert len(str1) == len(str2)
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

