"""The implementations under comparison.

A runner turns a built case directory into ``mageckGSEA`` output text. There
are two: the compiled helper, which is the oracle while it still exists, and
the Python port, which does not exist yet.

``PythonRunner`` is deliberately written before the code it calls, so that the
port's interface is decided by what a caller needs rather than inherited from
the shape of a command line. ``pathwayFunc.py:291-300`` builds a shell string,
runs it, and reads a file back; the point of the port is that
``mageck pathway`` should be able to call a function and get objects, the same
way the RRA port did.
"""

import hashlib
import os
import shutil
import subprocess


class RunnerUnavailable(Exception):
    """Raised when an implementation is not present to be tested."""


class CppRunner:
    """The compiled helper, invoked as ``pathwayFunc.py`` invokes it.

    Invoked by bare name so it picks up whichever build is first on PATH,
    matching how the package actually calls it.
    """

    name = "cpp"

    def __init__(self, binary="mageckGSEA"):
        self.binary = binary

    def available(self):
        return shutil.which(self.binary) is not None

    def fingerprint(self):
        """SHA-256 of the binary, so a golden can be traced to a build.

        Only the basename is recorded, not the resolved path: the hash is what
        identifies the build, while the absolute path would additionally bake
        the recording machine's home directory into a file committed to a
        public repository.
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
               "-g", str(directory / "sets.gmt"),
               "-r", str(directory / "rank.txt"),
               "-o", str(out),
               "-c", str(case.score_column),
               "-p", str(case.perms)]
        if case.sort_by_p:
            cmd.append("-s")
        if case.skip_header:
            cmd.append("-H")
        if case.reverse:
            cmd.append("-e")
        if case.pathway_name is not None:
            cmd += ["-n", case.pathway_name]
        if case.seed is not None:
            cmd += ["--seed", str(case.seed)]
        return cmd, out

    def run(self, case, directory):
        """Run the helper and return ``(output_text, log_text)``.

        Raises on a non-zero exit so a crash cannot be mistaken for a result.
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
                "mageckGSEA exited 0 but wrote no output for case %r"
                % case.name)
        return out.read_text(), proc.stdout + proc.stderr


class PythonRunner:
    """The Python port. Not yet implemented.

    The port is expected to expose::

        mageck2.gsea.enrich(pathways, genes, scores, perms=1000, seed=123456,
                            ...) -> GSEAResult

    taking a mapping of gene sets and two parallel sequences rather than two
    file paths, and returning an object with one row per gene set, so that
    ``mageck pathway`` can call it in process. A ``enrich_from_file`` helper
    keeps the file-shaped entry point for the harness and for the CLI, and
    ``to_text()`` renders the result in the historical format -- which is what
    lets goldens recorded from the C++ remain usable as the oracle after the
    C++ is gone.
    """

    name = "python"

    def available(self):
        try:
            import mageck2.gsea  # noqa: F401
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
            from mageck2 import gsea
        except ImportError as exc:
            raise RunnerUnavailable(
                "mageck2.gsea does not exist yet; this is the port's target "
                "interface, see PythonRunner's docstring") from exc
        kwargs = dict(
            score_column=case.score_column,
            perms=case.perms,
            sort_by_p=case.sort_by_p,
            skip_header=case.skip_header,
            reverse=case.reverse,
            pathway_name=case.pathway_name,
        )
        if case.seed is not None:
            kwargs["seed"] = case.seed
        result = gsea.enrich_from_file(
            str(directory / "sets.gmt"),
            str(directory / "rank.txt"),
            **kwargs)
        return result.to_text(), ""
