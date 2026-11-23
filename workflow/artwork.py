"""Export the actual embedded figure images from author-supplied clean Word files."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import posixpath
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile

NS = {"w":"http://schemas.openxmlformats.org/wordprocessingml/2006/main",
      "a":"http://schemas.openxmlformats.org/drawingml/2006/main",
      "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def embedded(doc, prefix, total):
    selected={}
    with ZipFile(doc) as archive:
        rels=ET.fromstring(archive.read("word/_rels/document.xml.rels"))
        targets={r.get("Id"):r.get("Target") for r in rels if r.get("TargetMode")!="External"}
        tree=ET.fromstring(archive.read("word/document.xml"))
        previous=[]
        for p in tree.findall(".//w:body/w:p",NS):
            text="".join(t.text or "" for t in p.findall(".//w:t",NS)).strip()
            images=[item.get("{"+NS["r"]+"}embed") for item in p.findall(".//a:blip",NS)]
            caption=re.match(r"^Figure\s+("+prefix+r"\d+)\s*[.:]",text)
            if caption:
                fig=caption.group(1)
                matches=images or previous
                if len(matches)!=1 or fig in selected:
                    raise ValueError(f"Cannot assign a unique image to Figure {fig}")
                target=posixpath.normpath(posixpath.join("word",targets[matches[0]]))
                if not target.startswith("word/media/"):
                    raise ValueError(f"Unexpected embedded media relationship for Figure {fig}")
                selected[fig]=(PurePosixPath(target).suffix,archive.read(target))
                previous=[]
            elif images:
                previous=images
            elif text:
                previous=[]
    expected={prefix+str(i) for i in range(1,total+1)}
    if set(selected)!=expected:raise ValueError(f"Unexpected figure set in {doc.name}: {set(selected)^expected}")
    return selected


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("main","si","out"):p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    panels={**embedded(a.main,"",5),**embedded(a.si,"S",28)}
    a.out.mkdir(parents=True,exist_ok=False)
    hashes={}
    for fig,(suffix,data) in panels.items():
        name=f"Figure_{fig}{suffix}"
        (a.out/name).write_bytes(data)
        hashes[name]=hashlib.sha256(data).hexdigest()
    (a.out/"checksums.json").write_text(json.dumps(hashes,indent=2))
    print(f"Exported {len(panels)} caption-matched original images")


if __name__=="__main__":
    main()
