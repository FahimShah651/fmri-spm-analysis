import urllib.request, sys, os
urls_file = sys.argv[1]
dest_dir = sys.argv[2]
import json
with open(urls_file) as f:
    tasks = json.load(f)
results = {}
for task in tasks:
    dest = os.path.join(dest_dir, task["file"])
    success = False
    for url in task["urls"]:
        try:
            urllib.request.urlretrieve(url, dest)
            sz = os.path.getsize(dest)
            if sz > 10000:
                results[task["file"]] = {"status":"OK","size":sz}
                success = True
                break
        except:
            pass
    if not success:
        results[task["file"]] = {"status":"FAIL","size":0}
print(json.dumps(results))
