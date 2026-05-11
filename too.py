import time
import threading
import tkinter as tk
from tkinter import ttk, messagebox
from collections import defaultdict

try:
    import psutil
except ImportError:
    psutil = None

try:
    import pydivert
except ImportError:
    pydivert = None


def format_speed(bytes_per_second: float) -> str:
    if bytes_per_second < 1024:
        return f"{bytes_per_second:.0f} B/s"

    kb = bytes_per_second / 1024
    if kb < 1024:
        return f"{kb:.2f} KB/s"

    mb = kb / 1024
    if mb < 1024:
        return f"{mb:.2f} MB/s"

    gb = mb / 1024
    return f"{gb:.2f} GB/s"


def format_mbps(bytes_per_second: float) -> str:
    return f"{bytes_per_second * 8 / 1024 / 1024:.2f} Mbps"


class ProcessNetworkMonitor:
    def __init__(self):
        self.lock = threading.Lock()

        self.pid_bytes = defaultdict(lambda: {"down": 0, "up": 0})
        self.total_down = 0
        self.total_up = 0

        self.port_pid_map = {}
        self.running = True

        self.last_connection_refresh = 0

    def refresh_connections(self):
        """
        建立 本地端口 -> PID 的映射。
        抓到数据包后，根据端口反查是哪个进程。
        """
        mapping = {}

        try:
            connections = psutil.net_connections(kind="inet")
        except Exception:
            return

        for conn in connections:
            if not conn.pid:
                continue

            if not conn.laddr:
                continue

            try:
                local_ip = conn.laddr.ip
                local_port = conn.laddr.port
            except Exception:
                continue

            proto = "tcp" if conn.type == 1 else "udp"

            mapping[(local_ip, local_port, proto)] = conn.pid
            mapping[("0.0.0.0", local_port, proto)] = conn.pid
            mapping[("::", local_port, proto)] = conn.pid
            mapping[(local_port, proto)] = conn.pid

        with self.lock:
            self.port_pid_map = mapping

    def get_process_name(self, pid):
        try:
            return psutil.Process(pid).name()
        except Exception:
            return f"PID {pid}"

    def find_pid_for_packet(self, packet):
        proto = None

        if packet.tcp:
            proto = "tcp"
        elif packet.udp:
            proto = "udp"
        else:
            return None

        if packet.is_outbound:
            local_ip = packet.src_addr
            local_port = packet.src_port
        else:
            local_ip = packet.dst_addr
            local_port = packet.dst_port

        with self.lock:
            mapping = self.port_pid_map.copy()

        pid = mapping.get((local_ip, local_port, proto))
        if pid is not None:
            return pid

        pid = mapping.get(("0.0.0.0", local_port, proto))
        if pid is not None:
            return pid

        pid = mapping.get(("::", local_port, proto))
        if pid is not None:
            return pid

        return mapping.get((local_port, proto))

    def capture_loop(self):
        self.refresh_connections()

        try:
            with pydivert.WinDivert("tcp or udp") as w:
                for packet in w:
                    if not self.running:
                        break

                    now = time.time()
                    if now - self.last_connection_refresh > 1:
                        self.refresh_connections()
                        self.last_connection_refresh = now

                    pid = self.find_pid_for_packet(packet)
                    packet_size = len(packet.raw)

                    with self.lock:
                        if packet.is_outbound:
                            self.total_up += packet_size
                            if pid:
                                self.pid_bytes[pid]["up"] += packet_size
                        else:
                            self.total_down += packet_size
                            if pid:
                                self.pid_bytes[pid]["down"] += packet_size

                    w.send(packet)

        except Exception as e:
            print("Capture error:", e)

    def start(self):
        t = threading.Thread(target=self.capture_loop, daemon=True)
        t.start()

    def snapshot(self):
        with self.lock:
            return {
                "total_down": self.total_down,
                "total_up": self.total_up,
                "pid_bytes": {
                    pid: data.copy()
                    for pid, data in self.pid_bytes.items()
                }
            }


class NetSpeedWindow:
    def __init__(self, root):
        self.root = root
        self.root.title("实时带宽监控")
        self.root.geometry("620x360")
        self.root.minsize(560, 320)
        self.root.resizable(True, True)

        self.monitor = ProcessNetworkMonitor()
        self.monitor.start()

        self.last_snapshot = self.monitor.snapshot()
        self.last_time = time.time()

        self.build_ui()
        self.update_speed()

        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

    def build_ui(self):
        self.root.configure(padx=18, pady=14)

        title = ttk.Label(
            self.root,
            text="实时带宽监控",
            font=("Microsoft YaHei UI", 18, "bold")
        )
        title.pack(pady=(0, 12))

        self.total_label = ttk.Label(
            self.root,
            text="总下载：0.00 KB/s    0.00 Mbps    |    总上传：0.00 KB/s    0.00 Mbps",
            font=("Microsoft YaHei UI", 11)
        )
        self.total_label.pack(anchor="w", pady=(0, 10))

        columns = ("rank", "process", "download", "upload", "total")

        self.tree = ttk.Treeview(
            self.root,
            columns=columns,
            show="headings",
            height=6
        )

        self.tree.heading("rank", text="排名")
        self.tree.heading("process", text="进程 / 应用")
        self.tree.heading("download", text="下载")
        self.tree.heading("upload", text="上传")
        self.tree.heading("total", text="总占用")

        self.tree.column("rank", width=60, anchor="center")
        self.tree.column("process", width=180, anchor="w")
        self.tree.column("download", width=120, anchor="center")
        self.tree.column("upload", width=120, anchor="center")
        self.tree.column("total", width=120, anchor="center")

        self.tree.pack(fill="both", expand=True)

        note = ttk.Label(
            self.root,
            text="说明：此窗口显示当前电脑实时流量 Top 3。若没有数据，请用管理员身份运行 PyCharm。",
            font=("Microsoft YaHei UI", 9)
        )
        note.pack(anchor="w", pady=(10, 0))

    def update_speed(self):
        current_snapshot = self.monitor.snapshot()
        current_time = time.time()

        elapsed = current_time - self.last_time
        if elapsed <= 0:
            elapsed = 1

        total_down_speed = (
            current_snapshot["total_down"] - self.last_snapshot["total_down"]
        ) / elapsed

        total_up_speed = (
            current_snapshot["total_up"] - self.last_snapshot["total_up"]
        ) / elapsed

        self.total_label.config(
            text=(
                f"总下载：{format_speed(total_down_speed)}    {format_mbps(total_down_speed)}"
                f"    |    "
                f"总上传：{format_speed(total_up_speed)}    {format_mbps(total_up_speed)}"
            )
        )

        rows = []

        old_pid_bytes = self.last_snapshot["pid_bytes"]
        new_pid_bytes = current_snapshot["pid_bytes"]

        for pid, new_data in new_pid_bytes.items():
            old_data = old_pid_bytes.get(pid, {"down": 0, "up": 0})

            down_speed = max(0, (new_data["down"] - old_data["down"]) / elapsed)
            up_speed = max(0, (new_data["up"] - old_data["up"]) / elapsed)
            total_speed = down_speed + up_speed

            if total_speed <= 0:
                continue

            process_name = self.monitor.get_process_name(pid)

            rows.append({
                "pid": pid,
                "name": process_name,
                "down": down_speed,
                "up": up_speed,
                "total": total_speed
            })

        rows.sort(key=lambda x: x["total"], reverse=True)
        top_rows = rows[:3]

        for item in self.tree.get_children():
            self.tree.delete(item)

        for index, row in enumerate(top_rows, start=1):
            display_name = f"{row['name']}  ({row['pid']})"

            self.tree.insert(
                "",
                "end",
                values=(
                    index,
                    display_name,
                    format_speed(row["down"]),
                    format_speed(row["up"]),
                    format_speed(row["total"])
                )
            )

        self.last_snapshot = current_snapshot
        self.last_time = current_time

        self.root.after(1000, self.update_speed)

    def on_close(self):
        self.monitor.running = False
        self.root.destroy()


def main():
    if psutil is None:
        root = tk.Tk()
        root.withdraw()
        messagebox.showerror(
            "缺少依赖",
            "缺少 psutil。\n\n请在 PyCharm 的 Python Packages 里搜索 psutil 并安装。"
        )
        return

    if pydivert is None:
        root = tk.Tk()
        root.withdraw()
        messagebox.showerror(
            "缺少依赖",
            "缺少 pydivert。\n\n请在 PyCharm 的 Python Packages 里搜索 pydivert 并安装。"
        )
        return

    root = tk.Tk()
    NetSpeedWindow(root)
    root.mainloop()


if __name__ == "__main__":
    main()
