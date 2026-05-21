"""
Auto-build a PowerPoint from ROOT/plot PDFs in EOS.

Example:
    python3 make_muon_jet_slides.py \
        --out /eos/user/d/dblend/PowerPoint/Muons_Jets_KIT_noScaleSmear_noJERC_noJER_BTagRW_ZPtRW.pptx

Optional:
  --layout single   # one PDF per slide
  --layout grid     # 2x2 PDFs per slide, default
  --recursive       # include PDFs in subdirectories too
"""
#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import tempfile
import os
from pathlib import Path
from typing import List, Tuple

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor

ERAS = ["12022", "22022", "12023", "22023", "2024", "2025"]
PLOT_GROUPS = [

    # Standard overlay PDFs
    {
        "name": "Muons",
        "type": "pdf",
        "path": "Muons/KIT/noScaleSmear_noJERC_noJER_BTagRW_ZPtRW",
    },

    {
        "name": "Jets",
        "type": "pdf",
        "path": "Jets/KIT/noScaleSmear_noJERC_noJER_BTagRW_ZPtRW",
    },

    # BTag efficiency PDFs
    {
        "name": "BTagEfficiency",
        "type": "pdf",
        "path": "BTagEfficiency/KIT/noScaleSmear_noJERC_noJER_noBTagRW_noZPtRW",
    },

    # Mass resolution single PNG
    {
        "name": "MassResolution",
        "type": "png",
        "path": "Mass_Resolutions/ZCR/KIT/noScaleSmear_noJERC_noJER_BTagRW_ZPtRW",
    },
]
CONFIG = "noScaleSmear_noJERC_noJER_BTagRW_ZPtRW"
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
        default="/eos/user/d/dblend/PowerPoint/Muons_Jets_KIT_noScaleSmear_noJERC_noJER_BTagRW_ZPtRW.pptx"
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
    subtitle.text = "KIT / BTagRW / ZPtRW"

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
