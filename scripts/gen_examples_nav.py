import pathlib

root = pathlib.Path(__file__).parent.parent
docs = root / "docs"
examples = docs / "examples"
notebooks = docs / "notebooks"

base = ["  - Examples:",]


def add_tree(base, folder, ext="py"):
    for file in sorted(folder.rglob(f"*.{ext}")):
        example_path = file.relative_to(docs)

        if file.parent.name.startswith("."):
            continue
        if file.name.startswith("_wip"):
            continue

        base.append(f"    - {file.stem}: {str(example_path)}")


add_tree(base, examples, ext="py")
base.append("  - Notebooks:")
add_tree(base, notebooks, ext="ipynb")

print("\n".join(base))
