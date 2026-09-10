"""
Auto-build a PowerPoint from ROOT/plot PDFs in EOS.

Example:
    python3 make_muon_jet_slides.py \
        --out /eos/user/d/dblend/PowerPoint/Muons_Jets_KIT_noScaleSmear_noBTagRW_noZPtRW.pptx

Optional:
  --layout single   # one PDF per slide
  --layout grid     # 2x2 PDFs per slide, default
  --recursive       # include PDFs in subdirectories too
"""
#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import tempfile
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor

ERAS = ["12022", "22022", "12023", "22023", "2024", "2025"]
ERA_LABELS = {
    "12022": "2022",
    "22022": "2022EE",
    "12023": "2023",
    "22023": "2023BPix",
    "2024": "2024",
    "2025": "2025",
}

SYST_ORDER = [
    "QCD",
    "PDF",
    "ISR",
    "FSR",
    "Pileup",
    "MuonSF",
    "Trigger",
    "JEC",
    "JER",
    "BTagSF",
    "BTagEffStat",
    "ZPt",
]

PLOT_GROUPS = [

    # Standard overlay PDFs
    {
        "name": "Muons",
        "type": "pdf",
        "path": "Muons/KIT/noScaleSmear_BTagRW_ZPtRW",
    },

    {
        "name": "Jets",
        "type": "pdf",
        "path": "Jets/KIT/noScaleSmear_BTagRW_ZPtRW",
    },

    # BTag efficiency PDFs
    {
        "name": "BTagEfficiency",
        "type": "pdf",
        "path": "BTagEfficiency/KIT/noScaleSmear_BTagRW_ZPtRW",
    },
    
    # Pileup reweighting PDFs
    {
        "name": "PileupReweighting",
        "type": "pdf",
        "path": "PileupReweighting/KIT/noScaleSmear_BTagRW_ZPtRW",
    },
    
    # Mass resolution single PNG
    {
        "name": "MassResolution",
        "type": "png",
        "path": "Mass_Resolutions/ZCR/KIT/noScaleSmear_BTagRW_ZPtRW",
    },
]
CONFIG = "noScaleSmear_BTagRW_ZPtRW"
CORR = "KIT"


def find_pdfs(folder: Path):
    return sorted(folder.glob("*.pdf"))


def pdf_to_png(pdf: Path, outdir: Path, dpi: int = 180) -> Path:

    out_png = outdir / f"{pdf.stem}.png"

    import fitz

    doc = fitz.open(str(pdf))
    page = doc.load_page(0)

    zoom = dpi / 72.0
    mat = fitz.Matrix(zoom, zoom)

    pix = page.get_pixmap(matrix=mat, alpha=False)
    pix.save(str(out_png))

    doc.close()

    return out_png

def add_text(slide, text, x, y, w, h, size=12, bold=False,
             color=None, align=PP_ALIGN.LEFT):

    box = slide.shapes.add_textbox(
        Inches(x), Inches(y), Inches(w), Inches(h)
    )

    tf = box.text_frame
    tf.clear()

    p = tf.paragraphs[0]
    p.text = text
    p.alignment = align

    run = p.runs[0]
    run.font.size = Pt(size)
    run.font.bold = bold

    if color is not None:
        run.font.color.rgb = color

    return box

def add_image(slide, img, x, y, w, h):

    slide.shapes.add_picture(
        str(img),
        Inches(x),
        Inches(y),
        width=Inches(w),
        height=Inches(h),
    )

def add_grid_slides(prs, plots, era, category):

    # Layout 14 = 3x2 Grid Layout
    layout = prs.slide_layouts[14]

    # picture placeholders
    picture_ids = [11, 12, 13, 14, 15, 16]

    for i in range(0, len(plots), 6):

        chunk = plots[i:i+6]

        slide = prs.slides.add_slide(layout)

        # title
        slide.shapes.title.text = f"{category} — {era}"

        for (img, pdf), pid in zip(chunk, picture_ids):

            try:
                placeholder = slide.placeholders[pid]

                slide.shapes.add_picture(
                    str(img),
                    placeholder.left,
                    placeholder.top,
                    width=placeholder.width,
                    height=placeholder.height,
                )

                # remove the empty placeholder box
                sp = placeholder._element
                sp.getparent().remove(sp)

            except Exception as e:
                print(f"Could not insert {img}")
                print(e)

def add_picture_fit_emu(slide, img, x, y, w, h):

    from PIL import Image

    with Image.open(img) as im:
        iw, ih = im.size

    box_ratio = w / h
    img_ratio = iw / ih

    if img_ratio > box_ratio:
        new_w = w
        new_h = int(w / img_ratio)
    else:
        new_h = h
        new_w = int(h * img_ratio)

    new_x = x + int((w - new_w) / 2)
    new_y = y + int((h - new_h) / 2)

    slide.shapes.add_picture(
        str(img),
        new_x,
        new_y,
        width=new_w,
        height=new_h,
    )



def pct_text(value, digits=2):
    """Format a percent value stored as a string/float from the CSV."""
    if value is None or value == "":
        return "--"
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def read_systematics_csv(csv_path: Path):
    """
    Read the CSV written by PrintSystematics_withCSV.C.

    If the CSV has multiple entries for the same era/source because the macro was
    re-run in append mode, keep the last row for each era/source.
    """
    single: Dict[Tuple[str, str], dict] = {}
    ranges: Dict[Tuple[str, str], dict] = {}

    if not csv_path.exists():
        print(f"[WARNING] Systematics CSV not found: {csv_path}")
        return single, ranges

    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            era = row.get("era", "").strip()
            source = row.get("source", "").strip()
            row_type = row.get("row_type", "").strip()
            if not era or not source:
                continue
            if row_type == "single":
                single[(era, source)] = row
            elif row_type == "range":
                ranges[(era, source)] = row

    return single, ranges


def add_table_slide(prs, title, headers, rows, footer=None, font_size=11):
    """Add a clean, PowerPoint-native table slide."""
    layout = prs.slide_layouts[7]  # Title Only
    slide = prs.slides.add_slide(layout)
    slide.shapes.title.text = title

    x = Inches(0.55)
    y = Inches(1.25)
    w = Inches(12.25)
    h = Inches(4.95)

    table_shape = slide.shapes.add_table(
        len(rows) + 1,
        len(headers),
        x,
        y,
        w,
        h,
    )
    table = table_shape.table

    # Column widths
    table.columns[0].width = Inches(1.65)
    if len(headers) > 1:
        other_w = int((w - Inches(1.65)) / (len(headers) - 1))
        for i in range(1, len(headers)):
            table.columns[i].width = other_w

    # Header
    for j, text in enumerate(headers):
        cell = table.cell(0, j)
        cell.text = text
        cell.fill.solid()
        cell.fill.fore_color.rgb = RGBColor(0, 0, 0)
        for p in cell.text_frame.paragraphs:
            p.alignment = PP_ALIGN.CENTER
            for run in p.runs:
                run.font.bold = True
                run.font.size = Pt(font_size)
                run.font.color.rgb = RGBColor(255, 255, 255)

    # Body
    for i, row in enumerate(rows, start=1):
        for j, text in enumerate(row):
            cell = table.cell(i, j)
            cell.text = str(text)
            for p in cell.text_frame.paragraphs:
                p.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
                for run in p.runs:
                    run.font.size = Pt(font_size)
                    run.font.bold = (j == 0)

    if footer:
        add_text(
            slide,
            footer,
            0.65,
            6.55,
            12.0,
            0.35,
            size=9,
            color=RGBColor(80, 80, 80),
            align=PP_ALIGN.LEFT,
        )

    return slide


def add_systematics_summary_slides(prs, csv_path: Optional[Path], include_details=True):
    if csv_path is None:
        return

    single, ranges = read_systematics_csv(csv_path)
    if not single and not ranges:
        return

    era_headers = [ERA_LABELS.get(e, e) for e in ERAS]

    # 1) Envelope comparison table
    headers = ["Source"] + era_headers
    rows = []
    for source in SYST_ORDER:
        row = [source]
        for era in ERAS:
            row.append(pct_text(single.get((era, source), {}).get("envelope_pct")))
        rows.append(row)

    add_table_slide(
        prs,
        "Systematic uncertainties: DiMuonMass_Higgs envelope (%)",
        headers,
        rows,
        footer=f"Read from {csv_path.name}; values are envelope = max(|up|, |down|).",
        font_size=9,
    )

    # 2) BDT input range comparison table
    range_rows = []
    for source in SYST_ORDER:
        row = [source]
        for era in ERAS:
            val = ranges.get((era, source), {}).get("range_text", "--")
            row.append(val if val else "--")
        range_rows.append(row)

    add_table_slide(
        prs,
        "Systematic uncertainty ranges across BDT inputs (%)",
        headers,
        range_rows,
        footer="Range is min–max envelope across the BDT variables listed in PrintSystematics.",
        font_size=10,
    )

    # 3) Optional one detailed table per era
    if include_details:
        for era in ERAS:
            detail_rows = []
            for source in SYST_ORDER:
                r = single.get((era, source), {})
                rr = ranges.get((era, source), {})
                detail_rows.append([
                    source,
                    pct_text(r.get("up_pct")),
                    pct_text(r.get("down_pct")),
                    pct_text(r.get("envelope_pct")),
                    rr.get("range_text", "--") or "--",
                ])

            add_table_slide(
                prs,
                f"Systematics detail — {ERA_LABELS.get(era, era)}",
                ["Source", "Up (%)", "Down (%)", "Envelope (%)", "BDT range (%)"],
                detail_rows,
                footer="Single histogram: DiMuonMass_Higgs in corr_DM_VBF_SR.",
                font_size=10,
            )

def add_era_comparison_slides(prs, tmpdir, base):

    layout = prs.slide_layouts[7]  # Title Only

    slide_w = prs.slide_width
    slide_h = prs.slide_height

    left = Inches(0.65)
    right = Inches(0.65)
    top = Inches(1.55)
    bottom = Inches(1.05)

    h_gap = Inches(0.35)
    v_gap = Inches(0.45)

    usable_w = slide_w - left - right
    usable_h = slide_h - top - bottom

    box_w = (usable_w - 2*h_gap) / 3
    box_h = (usable_h - v_gap) / 2

    positions = []
    for row in range(2):
        for col in range(3):
            x = left + col * (box_w + h_gap)
            y = top + row * (box_h + v_gap)
            positions.append((x, y, box_w, box_h))

    for group in PLOT_GROUPS:

        category = group["name"]

        hist_names = set()

        for era in ERAS:
            folder = Path(base) / era / group["path"]
            if folder.exists():
                if group["type"] == "pdf":

                    for pdf in find_pdfs(folder):
                        hist_names.add(pdf.name)

                elif group["type"] == "png":
                    
                    for png in sorted(folder.glob("*.png")):
                        hist_names.add(png.name)

        for hist_name in sorted(hist_names):

            slide = prs.slides.add_slide(layout)
            slide.shapes.title.text = f"{category}: {hist_name.replace('.pdf','')}"

            # Footer text
            add_text(
                slide,
                "University of Iowa",
                3.15,
                7.15,
                2.5,
                0.25,
                size=9,
                bold=True,
                color=RGBColor(255, 255, 255),
            )

            add_text(
                slide,
                "May 21",
                6.35,
                7.15,
                1.2,
                0.25,
                size=9,
                bold=True,
                color=RGBColor(255, 255, 255),
                align=PP_ALIGN.CENTER,
            )

            for era, (x, y, w, h) in zip(ERAS, positions):

                plot_file = Path(base) / era / group["path"] / hist_name

                if not plot_file.exists():
                    continue

                local_tmp = tmpdir / "era_compare" / category / hist_name.replace(".pdf", "")
                local_tmp.mkdir(parents=True, exist_ok=True)

                try:
                    if group["type"] == "pdf":
                        img = pdf_to_png(plot_file, local_tmp)

                    else:
                        img = plot_file

                    # outlined era label above plot
                    label_box = slide.shapes.add_shape(
                        1,
                        x,
                        y - Inches(0.22),
                        Inches(0.65),
                        Inches(0.16),
                    )

                    # white fill
                    label_box.fill.solid()
                    label_box.fill.fore_color.rgb = RGBColor(255, 255, 255)

                    # yellow outline
                    label_box.line.color.rgb = RGBColor(255, 205, 0)
                    label_box.line.width = Pt(1.5)

                    add_text(
                        slide,
                        era,
                        x / 914400.0 + 0.015,
                        y / 914400.0 - 0.215,
                        0.58,
                        0.14,
                        size=7,
                        bold=True,
                        color=RGBColor(0, 0, 0),
                        align=PP_ALIGN.CENTER,
                    )

                    add_picture_fit_emu(slide, img, x, y, w, h)

                except Exception as e:
                    print(f"Failed era comparison: {pdf}")
                    print(e)

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--base",
        default="/eos/user/d/dblend/StackedPlots"
    )

    parser.add_argument(
        "--out",
        default="/eos/user/d/dblend/PowerPoint/Muons_Jets_KIT_noScaleSmear_BTagRW_ZPtRW.pptx"
    )

    parser.add_argument(
        "--syst-csv",
        default="/eos/user/d/dblend/PowerPoint/systematics_summary.csv",
        help="CSV written by PrintSystematics_withCSV.C. If missing, table slides are skipped.",
    )

    parser.add_argument(
        "--no-syst-detail-slides",
        action="store_true",
        help="Only add the two all-era summary systematics slides; skip one detailed table slide per era.",
    )

    args = parser.parse_args()

    # Load UIowa template
    prs = Presentation("/eos/user/d/dblend/PowerPoint/IOWA-PPT-Template-Standard-2025.pptx")

    # Remove existing demo/template slides
    xml_slides = prs.slides._sldIdLst

    for i in range(len(prs.slides)-1, -1, -1):
        rId = xml_slides[i].rId
        prs.part.drop_rel(rId)
        del xml_slides[i]

    # ---------------------------------
    # Title slide
    # ---------------------------------
    title_slide = prs.slides.add_slide(prs.slide_layouts[0])

    title_slide.shapes.title.text = "Muon & Jet Validation"

    subtitle = title_slide.placeholders[1]
    subtitle.text = "KIT / no BTagRW / no ZPtRW" 

    # ---------------------------------
    # Systematics table slides
    # ---------------------------------
    syst_csv = Path(args.syst_csv) if args.syst_csv else None
    add_systematics_summary_slides(
        prs,
        syst_csv,
        include_details=not args.no_syst_detail_slides,
    )

    # ---------------------------------
    # Begin plots
    # ---------------------------------

    with tempfile.TemporaryDirectory() as tmp:

        tmpdir = Path(tmp)

        add_era_comparison_slides(
            prs,
            tmpdir,
            args.base,
        )

        local_out = Path("/tmp") / f"{Path(args.out).stem}_{os.getpid()}.pptx"

        prs.save(local_out)

        eos_out = Path(args.out)
        eos_out.parent.mkdir(parents=True, exist_ok=True)

        shutil.copy2(local_out, eos_out)

        print(f"[OK] Wrote local: {local_out}")
        print(f"[OK] Copied to EOS: {eos_out}")


if __name__ == "__main__":
    main()
