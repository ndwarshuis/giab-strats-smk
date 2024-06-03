from typing import TypeVar, Callable, NoReturn


X = TypeVar("X")
Y = TypeVar("Y")
Z = TypeVar("Z")


# 'id' in python apparently isn't this
def noop(x: X) -> X:
    return x


def compose(g: Callable[[Y], Z], f: Callable[[X], Y]) -> Callable[[X], Z]:
    return lambda x: g(f(x))


def from_maybe(default: X, x: X | None) -> X:
    return default if x is None else x


def fmap_maybe_def(default: Y, f: Callable[[X], Y], x: X | None) -> Y:
    return default if x is None else f(x)


def maybe_to_list(x: X | None) -> list[X]:
    return [] if x is None else [x]


def fmap_maybe(f: Callable[[X], Y], x: X | None) -> None | Y:
    return fmap_maybe_def(None, f, x)


def both(f: Callable[[X], Y], x: tuple[X, X]) -> tuple[Y, Y]:
    """2-tuple functor thingy"""
    return (f(x[0]), f(x[1]))


def thrice(f: Callable[[X], Y], x: tuple[X, X, X]) -> tuple[Y, Y, Y]:
    """3-tuple functor thingy"""
    return (f(x[0]), f(x[1]), f(x[2]))


def maybe2(x: tuple[X | None, X | None]) -> tuple[X, X] | None:
    if x[0] is not None and x[1] is not None:
        return (x[0], x[1])
    else:
        return None


def maybe3(x: tuple[X | None, X | None, X | None]) -> tuple[X, X, X] | None:
    if x[0] is not None and x[1] is not None and x[2] is not None:
        return (x[0], x[1], x[2])
    else:
        return None


def const(x: X) -> Callable[[Y], X]:
    def go(y: Y) -> X:
        return x

    return go


def constf(f: Callable[[X], Y]) -> Callable[[Z, X], Y]:
    return lambda _, x: f(x)


def seq(f: X, g: Y) -> Y:
    """The >> operator. For cases where I need more than on thing in a lambda."""
    f
    return g


def with_first(x: X, g: Callable[[X], Y]) -> X:
    """Kinda like finally but without the error handling magic.

    For cases where I want to compose to functions in a lambda but want the
    return value from the first.
    """
    g(x)
    return x


class DesignError(Exception):
    """Exception raised when the code is designed incorrectly (ie the 'this
    should not happen' error)"""

    pass


def raise_inline(msg: str) -> NoReturn:
    """Raise a design error.

    Useful for when I'm in a lambda and want to scream like a typical metalcore
    vocalist.
    """
    raise DesignError(msg)


def none_unsafe(x: X | None, f: Y, msg: None | str = None) -> Y:
    """Assert that something is None and continue if so."""
    if x is not None:
        raise DesignError(msg if msg is not None else "Should not be None")
    return f


def not_none_unsafe(x: X | None, f: Callable[[X], Y], msg: None | str = None) -> Y:
    """Call function with a value if it is not None, and error otherwise."""
    if x is None:
        raise DesignError(msg if msg is not None else "Should not be None")
    return f(x)


def match1_unsafe(xs: list[X], f: Callable[[X], Y], msg: None | str = None) -> Y:
    """Call function with the one value from a singleton list.

    Error if list is not a singleton.
    """
    match xs:
        case [x]:
            return f(x)
        case _:
            raise DesignError(
                msg if msg is not None else f"One input expected, got {len(xs)}"
            )


def match2_unsafe(xs: list[X], f: Callable[[X, X], Y], msg: None | str = None) -> Y:
    """Call function with the twos value from a 2-ary list.

    Error if list does not have two members.
    """
    match xs:
        case [x1, x2]:
            return f(x1, x2)
        case _:
            raise DesignError(
                msg if msg is not None else f"Two inputs expected, got {len(xs)}"
            )


def match12_unsafe(
    xs: list[X],
    f1: Callable[[X], Y],
    f2: Callable[[X, X], Y],
    msg: None | str = None,
) -> Y:
    """Combination of `match1_unsafe` and `match2_unsafe` with two functions
    for each case. Error if input list is does not have one or two elements.
    """
    match xs:
        case [x1]:
            return f1(x1)
        case [x1, x2]:
            return f2(x1, x2)
        case _:
            raise DesignError(
                msg if msg is not None else f"Two inputs expected, got {len(xs)}"
            )


# versions of unzip that won't return an empty tuple when given an empty list
def unzip2(xs: list[tuple[X, Y]]) -> tuple[list[X], list[Y]]:
    return ([x[0] for x in xs], [x[1] for x in xs])


def unzip3(xs: list[tuple[X, Y, Z]]) -> tuple[list[X], list[Y], list[Z]]:
    return ([x[0] for x in xs], [x[1] for x in xs], [x[2] for x in xs])


def filter_dict_strict(ds: dict[str, X], xs: set[str]) -> dict[str, X]:
    notindict = set(xs) - set(ds)
    if len(notindict) > 0:
        raise DesignError(f"Items not in input dictionary keys: {notindict}")
    return {k: v for k, v in ds.items() if k in xs}


def uncons_maybe(xs: list[X]) -> tuple[X, list[X]] | None:
    if len(xs) == 0:
        return None
    else:
        return (xs[0], xs[1:])
