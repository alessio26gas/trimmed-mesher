import subprocess
import os
import sys
import threading
from utils import resource_path


class ProcessHandler:
    def __init__(self, output_callback):
        self.process = None
        self.stop_reading = False
        self.output_callback = output_callback

    def start_process(self, inputs):
        if self.process and self.process.poll() is None:
            return

        exe = "mesher.exe" if sys.platform == "win32" else "mesher"
        command = [resource_path(exe)] + inputs   

        try:
            self.process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            threading.Thread(target=self._read_output, daemon=True).start()

        except Exception as e:
            self.output_callback(f"An error has occurred: {e}\n")

    def _read_output(self):
        if not self.process or not self.process.stdout:
            return

        for line in iter(self.process.stdout.readline, ''):
            if self.stop_reading:
                break
            self.output_callback(line)

    def stop_process(self):
        if self.process and self.process.poll() is None:
            self.output_callback("Aborting...\n")
            self.process.terminate()
            self.process.wait()
            self.output_callback("Trimmed Mesher has been interrupted.\n")
        else:
            self.output_callback("Trimmed Mesher is not currently running.\n")

    def cleanup(self):
        self.stop_reading = True
        if self.process and self.process.poll() is None:
            self.process.terminate()
            self.process.wait()