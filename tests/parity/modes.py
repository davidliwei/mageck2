"""How a corpus case is checked.

Shared by both ports, because the distinction is about the *status of the C++
behaviour*, not about which helper is being replaced.
"""

import enum


class Mode(enum.Enum):
    """Whether a case is checked against the C++, a contract, or a known fix.

    ``PARITY``
        The C++ behaviour is well defined and the port must reproduce it. This
        is the default and covers almost every case.

    ``CONTRACT``
        The case reaches genuinely undefined behaviour in the C++ -- an
        out-of-bounds read, an unguarded division -- so there is no correct
        answer to copy. Recording a golden would freeze whatever the compiler
        happened to emit on the recording machine into a requirement. These
        cases instead assert the behaviour the port is *choosing*, and the
        choice has to be made deliberately and documented in the test.

    ``DIVERGENCE``
        The C++ behaviour is well defined, reproducible, and wrong, and the
        port is required to differ from it. This is a separate mode rather
        than a flavour of ``CONTRACT`` because the distinction matters in both
        directions: the behaviour *is* defined, so it can be recorded and the
        record is worth having -- it documents what users of past versions
        actually got -- but it must not become a requirement, because the whole
        point of the case is that the port fixes it.

        A golden is still recorded for a ``DIVERGENCE`` case. The test then
        asserts two things: that the C++ still behaves as recorded (so the
        premise of the fix has not silently changed), and that the port does
        not, in the specific way the case names.
    """

    PARITY = "parity"
    CONTRACT = "contract"
    DIVERGENCE = "divergence"
