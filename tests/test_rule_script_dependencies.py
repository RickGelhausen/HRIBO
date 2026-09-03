"""Static guarantees for repository code executed by Snakemake rules."""

from __future__ import annotations

import ast
import re
from collections.abc import Iterator
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
RULES = REPO / "workflow" / "rules"
SCRIPTS = REPO / "workflow" / "scripts"
CODE_SUFFIXES = {".py", ".R", ".sh"}

RULE_START = re.compile(r"^rule\s+(?P<name>[A-Za-z_]\w*):", re.MULTILINE)
SCRIPT_PATH = re.compile(
    r"SCRIPTS(?P<parts>(?:\s*/\s*[\"'][^\"']+[\"'])+)"
)
NAMED_SCRIPT_INPUT = re.compile(
    r"(?P<key>[A-Za-z_]\w*)\s*=\s*str\(\s*SCRIPTS"
    r"(?P<parts>(?:\s*/\s*[\"'][^\"']+[\"'])+)\s*\)"
)
QUOTED_EXECUTABLE_INPUT = re.compile(
    r"(?:python(?:3)?|Rscript|bash|sh)\s+"
    r"\{input\.(?P<key>[A-Za-z_]\w*):q\}"
)
ANY_EXECUTABLE_INPUT = re.compile(
    r"(?:python(?:3)?|Rscript|bash|sh)\s+"
    r"\{input(?P<selector>\.[A-Za-z_]\w*|\[[^]]+\])?:q\}"
)


def _rule_bodies(text: str) -> Iterator[tuple[str, str]]:
    matches = list(RULE_START.finditer(text))
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        yield match.group("name"), text[match.start() : end]


def _directive_body(text: str, directive: str) -> str:
    """Return the indented body belonging to one Snakemake directive."""

    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip() != f"{directive}:":
            continue
        indentation = len(line) - len(line.lstrip())
        body = []
        for candidate in lines[index + 1 :]:
            if candidate.strip():
                candidate_indent = len(candidate) - len(candidate.lstrip())
                if candidate_indent <= indentation:
                    break
            body.append(candidate)
        return "\n".join(body)
    return ""


def _script_path(parts: str) -> Path:
    components = re.findall(r"/\s*[\"']([^\"']+)[\"']", parts)
    return SCRIPTS.joinpath(*components)


def _script_paths(text: str) -> set[Path]:
    return {_script_path(match.group("parts")) for match in SCRIPT_PATH.finditer(text)}


def _shared_script_inputs(text: str) -> dict[str, set[Path]]:
    """Find complete uppercase Python assignments that collect script paths."""

    shared = {}
    line_ends = [match.end() for match in re.finditer(r".*(?:\n|$)", text)]
    for match in re.finditer(r"^(?P<name>[A-Z][A-Z0-9_]*)\s*=", text, re.MULTILINE):
        # Rule files contain Snakemake directives that are not Python syntax.
        # Parse just the shortest complete assignment beginning at this match.
        for end in line_ends:
            if end <= match.start():
                continue
            candidate = text[match.start() : end]
            try:
                tree = ast.parse(candidate)
            except SyntaxError:
                continue
            if len(tree.body) != 1 or not isinstance(tree.body[0], ast.Assign):
                break
            paths = _script_paths(candidate)
            if paths:
                shared[match.group("name")] = paths
            break
    return shared


def _local_modules() -> dict[str, Path]:
    modules = {}
    for path in SCRIPTS.rglob("*.py"):
        relative = path.relative_to(SCRIPTS)
        module_parts = (
            relative.parent.parts
            if relative.name == "__init__.py"
            else relative.with_suffix("").parts
        )
        modules[".".join(module_parts)] = path
    return modules


def _resolved_module_paths(module: str, modules: dict[str, Path]) -> set[Path]:
    """Resolve a local module and the package initializers Python executes."""

    resolved = set()
    if module in modules:
        resolved.add(modules[module])
    parts = module.split(".")
    for length in range(1, len(parts)):
        package = ".".join(parts[:length])
        if package in modules:
            resolved.add(modules[package])
    return resolved


def _local_imports(path: Path, modules: dict[str, Path]) -> set[Path]:
    imports = set()
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.update(_resolved_module_paths(alias.name, modules))
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            imports.update(_resolved_module_paths(node.module, modules))
            # ``from lib import theme`` imports the local ``lib.theme`` module,
            # while ``from lib.io import read_genome`` only imports an attribute.
            for alias in node.names:
                imports.update(
                    _resolved_module_paths(f"{node.module}.{alias.name}", modules)
                )
    return imports


def _local_import_closure(path: Path, modules: dict[str, Path]) -> set[Path]:
    closure = set()
    pending = [path]
    while pending:
        for dependency in _local_imports(pending.pop(), modules):
            if dependency != path and dependency not in closure:
                closure.add(dependency)
                pending.append(dependency)
    return closure


def test_rule_executables_and_local_import_closures_are_job_inputs():
    """Keep repository code visible to Snakemake without leaking it into tools."""

    modules = _local_modules()
    executable_basenames = {
        path.name
        for path in SCRIPTS.rglob("*")
        if path.is_file() and path.suffix in CODE_SUFFIXES
    }
    failures = []
    inventory = []
    dependency_edges = set()

    for rule_path in sorted(RULES.glob("*.smk")):
        text = rule_path.read_text()
        shared_inputs = _shared_script_inputs(text)
        for rule_name, rule_body in _rule_bodies(text):
            input_body = _directive_body(rule_body, "input")
            shell_body = _directive_body(rule_body, "shell")
            if not shell_body:
                continue

            declared_code = _script_paths(input_body)
            for constant, paths in shared_inputs.items():
                if re.search(rf"\b{re.escape(constant)}\b", input_body):
                    declared_code.update(paths)

            named_executables = {
                match.group("key"): _script_path(match.group("parts"))
                for match in NAMED_SCRIPT_INPUT.finditer(input_body)
            }
            inventory.extend(
                (rule_path.name, rule_name, key, path)
                for key, path in named_executables.items()
            )

            for match in ANY_EXECUTABLE_INPUT.finditer(shell_body):
                selector = match.group("selector") or ""
                if not selector.startswith("."):
                    failures.append(
                        f"{rule_path.name}:{rule_name}: repository executable "
                        "input is positional rather than named"
                    )

            for key, executable in named_executables.items():
                quoted_invocations = {
                    match.group("key")
                    for match in QUOTED_EXECUTABLE_INPUT.finditer(shell_body)
                }
                if key not in quoted_invocations:
                    failures.append(
                        f"{rule_path.name}:{rule_name}: {executable.name} is not "
                        f"invoked through quoted input.{key}"
                    )
                if not executable.is_file():
                    failures.append(
                        f"{rule_path.name}:{rule_name}: missing executable {executable}"
                    )
                    continue
                if executable.suffix == ".py":
                    closure = _local_import_closure(executable, modules)
                    for importer in closure | {executable}:
                        dependency_edges.update(
                            (importer, imported)
                            for imported in _local_imports(importer, modules)
                        )
                    missing = sorted(closure - declared_code)
                    if missing:
                        failures.append(
                            f"{rule_path.name}:{rule_name}: undeclared local imports "
                            + ", ".join(
                                str(path.relative_to(SCRIPTS)) for path in missing
                            )
                        )

            if declared_code and "{input:q}" in shell_body:
                failures.append(
                    f"{rule_path.name}:{rule_name}: aggregate input passes code "
                    "dependencies to the command"
                )

            for basename in executable_basenames:
                literal = re.compile(
                    rf"(?<![A-Za-z0-9_.-]){re.escape(basename)}"
                    rf"(?![A-Za-z0-9_.-])"
                )
                if literal.search(shell_body):
                    failures.append(
                        f"{rule_path.name}:{rule_name}: literal repository "
                        f"executable {basename}"
                    )

    assert inventory, "no named repository executable inputs were discovered"
    assert (
        SCRIPTS / "excel_utils.py",
        SCRIPTS / "gff_utils.py",
    ) in dependency_edges
    assert any(
        dependency == SCRIPTS / "lib" / "__init__.py"
        for _, dependency in dependency_edges
    )
    assert failures == []
