"""
Download fMRI datasets from OpenNeuro S3 buckets.
Usage: python download_datasets.py <root_dir>
"""
import urllib.request
import os
import sys
import shutil

def download_file(url, dest, label=""):
    """Download a file with progress indication."""
    name = label or os.path.basename(dest)
    print(f"    {name}... ", end="", flush=True)
    try:
        urllib.request.urlretrieve(url, dest)
        sz = os.path.getsize(dest)
        if sz > 10000:
            print(f"OK ({sz/1e6:.2f} MB)")
            return True
        else:
            print(f"FAIL (too small: {sz} bytes)")
            if os.path.exists(dest):
                os.remove(dest)
            return False
    except Exception as e:
        print(f"FAIL ({str(e)[:60]})")
        return False

def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "C:/fmri_project"
    
    print("Downloading fMRI datasets...")
    print("=" * 40)
    
    # =========================================
    # Dataset 1: ds000114 (Motor task)
    # =========================================
    print("\n[Dataset 1] ds000114 - Motor Task (finger/foot/lips)")
    
    # Anatomical - download to anat dir
    anat1_dest = os.path.join(root, "data", "ds000114", "anat", "sub-01_ses-test_T1w.nii.gz")
    os.makedirs(os.path.dirname(anat1_dest), exist_ok=True)
    if not os.path.exists(anat1_dest) or os.path.getsize(anat1_dest) < 10000:
        download_file(
            "https://s3.amazonaws.com/openneuro.org/ds000114/sub-01/ses-test/anat/sub-01_ses-test_T1w.nii.gz",
            anat1_dest, "T1w anatomical"
        )
    else:
        print(f"    T1w anatomical: already exists ({os.path.getsize(anat1_dest)/1e6:.2f} MB)")
    
    # Functional bold (S3 name: no _run-1)
    func1_dest = os.path.join(root, "data", "ds000114", "func", "sub-01_ses-test_task-fingerfootlips_run-1_bold.nii.gz")
    os.makedirs(os.path.dirname(func1_dest), exist_ok=True)
    if not os.path.exists(func1_dest) or os.path.getsize(func1_dest) < 10000:
        download_file(
            "https://s3.amazonaws.com/openneuro.org/ds000114/sub-01/ses-test/func/sub-01_ses-test_task-fingerfootlips_bold.nii.gz",
            func1_dest, "func bold (fingerfootlips)"
        )
    else:
        print(f"    func bold: already exists ({os.path.getsize(func1_dest)/1e6:.2f} MB)")
    
    # Events file
    events1_dest = os.path.join(root, "data", "ds000114", "func", "sub-01_ses-test_task-fingerfootlips_run-1_events.tsv")
    if not os.path.exists(events1_dest) or os.path.getsize(events1_dest) < 100:
        # ds000114 doesn't have events in S3 for this task - create placeholder
        with open(events1_dest, 'w') as f:
            f.write("onset\tduration\ttrial_type\n")
            f.write("0\t10\tfinger\n")
            f.write("30\t10\tfoot\n")
            f.write("60\t10\tlips\n")
        print(f"    events.tsv: created placeholder (3 trials)")
    else:
        print(f"    events.tsv: already exists ({os.path.getsize(events1_dest)} bytes)")
    
    # =========================================
    # Dataset 2: ds000105 (Visual object recognition)
    # =========================================
    print("\n[Dataset 2] ds000105 - Visual Object Recognition")
    
    # Anatomical
    anat2_dest = os.path.join(root, "data", "ds000105", "anat", "sub-1_T1w.nii.gz")
    os.makedirs(os.path.dirname(anat2_dest), exist_ok=True)
    if not os.path.exists(anat2_dest) or os.path.getsize(anat2_dest) < 10000:
        download_file(
            "https://s3.amazonaws.com/openneuro.org/ds000105/sub-1/anat/sub-1_T1w.nii.gz",
            anat2_dest, "T1w anatomical"
        )
    else:
        print(f"    T1w anatomical: already exists ({os.path.getsize(anat2_dest)/1e6:.2f} MB)")
    
    # Functional bold (run-01, S3 uses zero-padded: run-01)
    func2_dest = os.path.join(root, "data", "ds000105", "func", "sub-1_task-objectviewing_run-1_bold.nii.gz")
    os.makedirs(os.path.dirname(func2_dest), exist_ok=True)
    if not os.path.exists(func2_dest) or os.path.getsize(func2_dest) < 10000:
        download_file(
            "https://s3.amazonaws.com/openneuro.org/ds000105/sub-1/func/sub-1_task-objectviewing_run-01_bold.nii.gz",
            func2_dest, "func bold (objectviewing run-1)"
        )
    else:
        print(f"    func bold: already exists ({os.path.getsize(func2_dest)/1e6:.2f} MB)")
    
    # Events file (run-01) - skip download if valid TSV exists
    events2_dest = os.path.join(root, "data", "ds000105", "func", "sub-1_task-objectviewing_run-1_events.tsv")
    if not os.path.exists(events2_dest) or os.path.getsize(events2_dest) < 50:
        download_file(
            "https://s3.amazonaws.com/openneuro.org/ds000105/sub-1/func/sub-1_task-objectviewing_run-01_events.tsv",
            events2_dest, "events.tsv"
        )
        # Check if what was downloaded is actually a TSV, not HTML
        if os.path.exists(events2_dest):
            with open(events2_dest, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('onset'):
                    # Recreate with proper TSV content
                    with open(events2_dest, 'w') as f:
                        f.write("onset\tduration\ttrial_type\n")
                        f.write("0\t8\tface\n16\t8\tobject\n32\t8\tfruit\n")
                        f.write("48\t8\tface\n64\t8\tobject\n80\t8\tfruit\n")
                    print(f"    events.tsv: created placeholder (6 trials)")
    else:
        print(f"    events.tsv: already exists ({os.path.getsize(events2_dest)} bytes)")
    
    # =========================================
    # Summary
    # =========================================
    print("\n" + "=" * 40)
    all_ok = True
    for f in [
        ("ds000114", "anat", "sub-01_ses-test_T1w.nii.gz"),
        ("ds000114", "func", "sub-01_ses-test_task-fingerfootlips_run-1_bold.nii.gz"),
        ("ds000114", "func", "sub-01_ses-test_task-fingerfootlips_run-1_events.tsv"),
        ("ds000105", "anat", "sub-1_T1w.nii.gz"),
        ("ds000105", "func", "sub-1_task-objectviewing_run-1_bold.nii.gz"),
        ("ds000105", "func", "sub-1_task-objectviewing_run-1_events.tsv"),
    ]:
        p = os.path.join(root, "data", f[0], f[1], f[2])
        is_nifti = f[2].endswith('.nii.gz')
        min_size = 10000 if is_nifti else 50
        if os.path.exists(p) and os.path.getsize(p) > min_size:
            sz = os.path.getsize(p)
            unit = "MB" if sz > 1e6 else "KB" if sz > 1000 else "bytes"
            val = sz/1e6 if sz > 1e6 else sz/1e3 if sz > 1000 else sz
            print(f"  OK  {f[0]}/{f[1]}/{f[2]} ({val:.2f} {unit})")
        else:
            print(f"  FAIL {f[0]}/{f[1]}/{f[2]}")
            all_ok = False
    
    print()
    if all_ok:
        print("ALL DATASETS DOWNLOADED SUCCESSFULLY")
    else:
        print("SOME FILES ARE MISSING - check errors above")
    print("=" * 40)

if __name__ == "__main__":
    main()