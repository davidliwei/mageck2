"""How a corpus case is checked.

Shared by both ports, because the distinction is about the *status of the C++
behaviour*, not about which helper is being replaced.
"""

import enum


class Mode(enum.Enum):
    """Whether a case is checked against the C++ or against a stated contract.

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

    There used to be a third mode, ``DIVERGENCE``, for behaviour that was
    well defined, reproducible, and wrong, where the port was required to
    differ. Five ``mageckGSEA`` input-handling bugs were classified that way.
    They were fixed in the C++ instead, which is the better answer when the
    C++ is still what ships: the cases became ordinary ``PARITY`` cases whose
    inputs happen to be malformed, and both implementations now have to agree
    with each other *and* with the repaired input. The mode is gone rather
    than kept empty because nothing needs it, and a mode no case uses is a
    claim the corpus does not support.
    """

    PARITY = "parity"
    CONTRACT = "contract"
