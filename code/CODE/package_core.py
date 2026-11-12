"""Small, tested author-review packaging utilities; no biological inference."""
import csv
import hashlib
import json
from pathlib import Path
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph
from openpyxl import Workbook, load_workbook


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    fields = fields or list(rows[0])
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def replace_paragraph(paragraph, expected, new):
    if paragraph.text != expected:
        raise ValueError('Paragraph source changed')
    paragraph.text = new


def insert_text(anchor, text, style='Normal'):
    element = OxmlElement('w:p')
    anchor._p.addprevious(element)
    p = Paragraph(element, anchor._parent)
    p.style = style
    p.add_run(text)
    return p


def safe_output(root, output, master):
    root, output, master = [Path(x).resolve() for x in (root, output, master)]
    if output == master or not output.is_relative_to(root):
        raise ValueError('Unsafe output target')
    return output


def classify_rates(rows):
    valid = [r for r in rows if r['Ka_Ks'] != 'NA']
    return {'computed': len(valid), 'withheld': len(rows)-len(valid),
            'scientific_status': 'REVIEW',
            'ratio_min': min((float(r['Ka_Ks']) for r in valid), default=None),
            'ratio_max': max((float(r['Ka_Ks']) for r in valid), default=None)}


def workbook(path, sheets):
    wb = Workbook()
    wb.remove(wb.active)
    expected = {}
    for name, rows in sheets.items():
        ws = wb.create_sheet(name[:31])
        fields = list(rows[0])
        values = [fields] + [[str(r.get(k, '')) for k in fields] for r in rows]
        for row in values:
            ws.append(row)
        ws.freeze_panes = 'A2'
        ws.auto_filter.ref = ws.dimensions
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(45, max(16, len(str(col[0].value))+2))
        expected[ws.title] = values
    wb.save(path)
    check = load_workbook(path, read_only=True, data_only=False)
    for name, values in expected.items():
        actual = [[str(v) if v is not None else '' for v in r] for r in check[name].values]
        assert actual == values, f'Workbook cell mismatch: {name}'
    check.close()


def manifest(directory):
    directory = Path(directory)
    entries = [(sha(p), str(p.relative_to(directory))) for p in sorted(directory.rglob('*'))
               if p.is_file() and p.name != 'SHA256SUMS.txt' and '__pycache__' not in p.parts]
    (directory / 'SHA256SUMS.txt').write_text(''.join(f'{digest}  {path}\n' for digest, path in entries))
    return len(entries)
