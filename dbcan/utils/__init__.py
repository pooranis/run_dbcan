import datetime
import sys

from .simplify_cgc import simplify_output
from .CGCFinder import cgc_finder



def printmsg(*args, **kwargs):
    """wrapper for print which prepends datetime and prints to STDERR, flushing immediately.

    Same args as :py:func:`print`, with addition of:
      begin (str): extra argument for adding string (like newlines) before printing line.

    """
    begin = kwargs.pop("begin", None)
    if begin is not None:
        print(begin, **kwargs, file=sys.stderr, flush=True)
    print("[", datetime.datetime.now(), "]", *args, **kwargs, file=sys.stderr, flush=True)

