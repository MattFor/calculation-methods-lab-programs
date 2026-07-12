from __future__ import annotations

import re
import sys
import csv
import math
import argparse
from pathlib import Path

try:
	from odf.opendocument import OpenDocumentSpreadsheet
	from odf.table import Table, TableRow, TableCell
	from odf.text import P
except ImportError as exc:
	raise SystemExit(
		"Missing dependency 'odfpy'. Install it with: pip install odfpy"
	) from exc

FORBIDDEN_SHEET_CHARS = re.compile(r"[\\/?*\[\]:]")


def sanitize_sheet_name(name: str, existing: set[str]) -> str:
	"""
	Return a Calc-safe sheet name, unique within the workbook.
	"""
	cleaned = FORBIDDEN_SHEET_CHARS.sub("_", name).strip()
	if not cleaned:
		cleaned = "Sheet"

	cleaned = cleaned[:31]
	base = cleaned
	suffix_i = 1
	while cleaned in existing:
		suffix = f"_{suffix_i}"
		cleaned = (base[: 31 - len(suffix)] + suffix)[:31]
		suffix_i += 1

	existing.add(cleaned)
	return cleaned


def parse_cell(value: str):
	"""
	Try to preserve numeric CSV cells as numbers; otherwise keep text.
	"""
	text = value.strip()
	if text == "":
		return "string", ""

	try:
		# Preserve integers and decimals for Calc formulas/charts.
		num = float(text)
		if math.isfinite(num):
			return "float", text
	except ValueError:
		pass

	return "string", value


def read_csv(path: Path) -> list[list[str]]:
	with path.open("r", newline="", encoding="utf-8-sig") as f:
		return [row for row in csv.reader(f)]


def add_csv_as_sheet(doc, sheet_name: str, rows: list[list[str]]) -> None:
	table = Table(name=sheet_name)

	for row in rows:
		tr = TableRow()
		for value in row:
			cell_type, cell_value = parse_cell(value)
			if cell_type == "float":
				cell = TableCell(valuetype="float", value=cell_value)
				cell.addElement(P(text=cell_value))
			else:
				cell = TableCell(valuetype="string")
				cell.addElement(P(text=cell_value))
			tr.addElement(cell)
		table.addElement(tr)

	doc.spreadsheet.addElement(table)


def collect_csv_files(paths: list[Path]) -> list[Path]:
	if paths:
		return [p for p in paths if p.suffix.lower() == ".csv"]

	data_dir = Path(__file__).resolve().parent.parent
	return sorted(data_dir.glob("*.csv"))


def main(argv: list[str]) -> int:
	parser = argparse.ArgumentParser(description="Create an ODS workbook from CSV files.")
	parser.add_argument(
		"csv_files",
		nargs="*",
		type=Path,
		help="CSV files to import. If omitted, uses data/*.csv",
	)

	parser.add_argument(
		"--output",
		type=Path,
		default=Path(__file__).resolve().parent.parent / "output" / "results.ods",
		help="Output ODS file (default: data/results.ods)",
	)

	args = parser.parse_args(argv)

	csv_files = collect_csv_files(args.csv_files)
	if not csv_files:
		raise SystemExit("No CSV files found.")

	doc = OpenDocumentSpreadsheet()
	existing_names: set[str] = set()

	for csv_path in csv_files:
		if not csv_path.exists():
			raise SystemExit(f"Missing file: {csv_path}")

		sheet_name = sanitize_sheet_name(csv_path.stem, existing_names)
		rows = read_csv(csv_path)
		add_csv_as_sheet(doc, sheet_name, rows)

	args.output.parent.mkdir(parents=True, exist_ok=True)
	doc.save(str(args.output))
	print(f"Saved: {args.output}")
	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv[1:]))
