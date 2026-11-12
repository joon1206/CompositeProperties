import csv

import pytest

from scripts.material_card import build_card


HEADER = "temperature_C,ef,nuf,af,vf,ep,nup,ap,phi,er,nur,ar,warp_crimp_deg,weft_crimp_deg,warp_share\n"
ROW = "25,82,0.22,2.49,0.42,70,0.17,0.5,0.70,2.5,0.37,50,5,15,0.5\n"


def test_card_writes_twelve_properties_per_input_temperature(tmp_path):
    source = tmp_path/"input.csv"
    target = tmp_path/"card.csv"
    source.write_text(HEADER+ROW+ROW.replace("25,", "100,", 1), encoding="utf-8")
    rows = build_card(source, target, "mori_tanaka")
    assert len(rows) == 2
    assert rows[0]["EX"] > rows[0]["EZ"]
    assert rows[0]["temperature_C"] == 25
    assert rows[1]["temperature_C"] == 100
    assert rows[0]["EX"] == pytest.approx(rows[1]["EX"])
    with target.open(newline="", encoding="utf-8") as f:
        saved = list(csv.DictReader(f))
    assert len(saved) == 2
    assert len([k for k in saved[0] if k in
                ("EX", "EY", "EZ", "GXY", "GYZ", "GXZ",
                 "NUXY", "NUYZ", "NUXZ", "AX", "AY", "AZ")]) == 12


def test_card_rejects_missing_input_instead_of_filling_it(tmp_path):
    source = tmp_path/"input.csv"
    source.write_text(HEADER.replace(",phi,", ",")+ROW, encoding="utf-8")
    with pytest.raises(ValueError, match="phi"):
        build_card(source, tmp_path/"card.csv", "legacy")


def test_card_accepts_direct_filled_matrix_measurements(tmp_path):
    source = tmp_path/"measured.csv"
    source.write_text(
        "temperature_C,ef,nuf,af,vf,warp_crimp_deg,weft_crimp_deg,warp_share,em_measured,num_measured,am_measured\n"
        "25,82,0.22,2.49,0.42,5,15,0.5,10,0.25,15\n", encoding="utf-8")
    rows = build_card(source, tmp_path/"card.csv", "measured")
    assert rows[0]["Em"] == pytest.approx(10)
    assert rows[0]["am"] == pytest.approx(15)
    assert rows[0]["EX"] > rows[0]["EZ"]
