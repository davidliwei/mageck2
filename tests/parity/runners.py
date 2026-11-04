"""The implementations under comparison.

A runner turns a built case directory into RRA output text. There are two: the
compiled helper, which is the oracle while it still exists, and the Python
port, which does not exist yet.

``PythonRunner`` is deliberately written now, before the code it calls. It
fixes the signature the port has to satisfy -- values in, a results table out,
no temporary files and no subprocess -- so the API shape is decided by the
harness rather than inherited from the shape of a command line. That is the
whole point of the port: ``crisprFunction.py:559-571`` builds a shell string
and reads a file back, which is why there is no in-process entry point today.
"""

import hashlib
import os
import shutil
import subprocess


class RunnerUnavailable(Exception):
    """Raised when an implementation is not present to be tested."""


class CppRunner:
    """The compiled RRA helper, invoked exactly as ``crisprFunction.py`` does.

    Invoked by bare name so it picks up whichever build is first on PATH,
    matching how the package actually calls it.
    """

    name = "cpp"

    def __init__(self, binary="RRA"):
        self.binary = binary

    def available(self):
        return shutil.which(self.binary) is not None

    def fingerprint(self):
        """SHA-256 of the binary, so a golden can be traced to a build.

        Only the basename is recorded, not the resolved path. The hash is what
        identifies the build; the absolute path would additionally bake the
        recording machine's home directory into a file committed to a public
        repository, which identifies nothing useful.
        """
        path = shutil.which(self.binary)
        if path is None:
            return None
        digest = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                digest.update(chunk)
        return "%s sha256:%s" % (os.path.basename(path), digest.hexdigest())

    def command(self, case, directory):
        out = directory / ("%s.out.txt" % self.name)
        cmd = [self.binary,
               "-i", str(directory / "input.txt"),
               "-o", str(out),
               "-p", str(case.alpha)]
        if case.uses_control:
            cmd += ["--control", str(directory / "control.txt")]
        cmd += list(case.flags)
        return cmd, out

    def run(self, case, directory):
        """Run the helper and return ``(output_text, stdout_text)``.

        Raises on a non-zero exit so a crash cannot be mistaken for a result,
        which matters for the ``CONTRACT`` cases where a crash is a plausible
        outcome.
        """
        if not self.available():
            raise RunnerUnavailable(
                "%r is not on PATH; run `pip install .` first" % self.binary)
        cmd, out = self.command(case, directory)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(
                proc.returncode, cmd, proc.stdout, proc.stderr)
        if not out.exists():
            raise RuntimeError(
                "RRA exited 0 but wrote no output for case %r" % case.name)
        return out.read_text(), proc.stdout + proc.stderr


class PythonRunner:
    """The Python port. Not yet implemented.

    The port is expected to expose::

        mageck2.rra.rank(items, groups, values, alpha=..., control=None,
                         permutation=100, by_group=True, ...) -> RRAResult

    taking arrays rather than a file path, and returning an object with one row
    per group, so that `mageck2 test` can call it in process. ``to_text()``
    renders that result in the historical file format, which is what lets the
    goldens recorded from the C++ remain usable as the oracle.
    """

    name = "python"

    def available(self):
        try:
            import mageck2.rra  # noqa: F401
        except ImportError:
            return False
        return True

    def fingerprint(self):
        try:
            import mageck2
        except ImportError:
            return None
        return "mageck2 %s" % getattr(mageck2, "__version__", "unknown")

    def run(self, case, directory):
        try:
            from mageck2 import rra
        except ImportError as exc:
            raise RunnerUnavailable(
                "mageck2.rra does not exist yet; this is the port's target "
                "interface, see PythonRunner's docstring") from exc
        control = None
        if case.uses_control:
            control = (directory / "control.txt").read_text().split()
        result = rra.rank_from_file(
            str(directory / "input.txt"),
            alpha=case.alpha,
            control=control,
            flags=list(case.flags),
        )
        return result.to_text(), ""
