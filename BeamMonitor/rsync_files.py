import subprocess
import sys

def rsync_single_file(remote_user, remote_host, remote_file, local_path):
    rsync_cmd = [
        "rsync",
        "-avz",
        "--progress",
        f"{remote_user}@{remote_host}:{remote_file}",
        local_path
    ]

    print(f"Syncing {remote_file}...")
    try:
        subprocess.run(rsync_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Rsync failed for {remote_file}: {e}")

if __name__ == "__main__":
    remote_user = "username"
    remote_host = "server.example.com"
    remote_base_folder = "/remote/folder/"
    local_path = "/mnt/ssd/backup/"

    txt_files = [
        "file1.txt",
        "file2.txt",
        "notes2025.txt"
    ]

    for filename in txt_files:
        remote_file = remote_base_folder + filename
        rsync_single_file(remote_user, remote_host, remote_file, local_path)

