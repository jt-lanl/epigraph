
class BestSoFar:
    '''keep track of the item with the highest value:
    Usage:
       so_far = BestSoFar()
       ...
       Loop:
          generate new item, compute its value
          so_far.update(value,item)
       ...
       so_far.best() is best item so far
       so_far.value is value of best item
    '''
    def __init__(self,value=None,item=None):
        self.item = item
        self.thevalue = value

    def _update_item(self,item):
        try:
            ## if it's copy-able, then make a copy
            self.item = item.copy()
        except AttributeError:
            self.item = item

    def update(self,value,item):
        if self.item is None or value > self.thevalue:
            self.thevalue = value
            self._update_item(item)
        return self

    @property
    def best(self):
        return self.item

    @property
    def value(self):
        return self.thevalue

