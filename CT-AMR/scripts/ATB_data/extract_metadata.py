import os
import time
import csv
import xml.etree.ElementTree as ET
from subprocess import run, PIPE

INPUT_FILE = "biosamples.txt"
OUTPUT_CSV = "1_biosample_metadata.csv"

def fetch_biosample_xml(biosample_id):
    cmd = f'esearch -db biosample -query "{biosample_id}" | efetch -format xml'
    result = run(cmd, shell=True, stdout=PIPE, stderr=PIPE, text=True)
    if result.returncode != 0 or not result.stdout:
        print(f"⚠️ Failed to fetch {biosample_id}")
        return None
    return result.stdout

def parse_xml(xml_content):
    metadata = {}
    try:
        root = ET.fromstring(xml_content)
        biosample = root.find(".//BioSample")
        if biosample is not None:
            metadata["BioSample"] = biosample.attrib.get("accession", "")
            for attr in biosample.findall(".//Attribute"):
                name = attr.attrib.get("attribute_name", "").strip()
                value = attr.text.strip() if attr.text else ""
                if name:
                    metadata[name] = value
    except ET.ParseError:
        print("⚠️ XML parse error")
    return metadata

def main():
    all_records = []
    with open(INPUT_FILE) as f:
        biosample_ids = [line.strip() for line in f if line.strip()]

    for biosample_id in biosample_ids:
        print(f"🔍 Fetching metadata for {biosample_id}...")
        xml = fetch_biosample_xml(biosample_id)
        if xml:
            record = parse_xml(xml)
            if record:
                all_records.append(record)
        time.sleep(0.4)  # Be polite to NCBI servers

    # Collect all unique field names
    all_keys = set()
    for r in all_records:
        all_keys.update(r.keys())
    all_keys = sorted(all_keys)

    # Write to CSV
    with open(OUTPUT_CSV, "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=all_keys)
        writer.writeheader()
        for row in all_records:
            writer.writerow(row)

    print(f"✅ Done! Metadata saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
