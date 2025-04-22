import itertools
import sys
import threading
import time

try:
    import emoji
except ImportError:
    emoji = None
import os

if not os.environ.get("PYTHONIOENCODING"):
    os.environ["PYTHONIOENCODING"] = "UTF-8"
    sys.stdout.reconfigure(encoding="utf-8")


class Spinner:
    """
    A colorful and animated loader for terminal applications.
    """

    def __init__(self, text="\nLoading...\n", spinner=None, color="cyan"):
        self.text = text
        self.spinner = spinner or ["⠋", "⠙", "⠸", "⠴", "⠦", "⠇"]
        self.spinner = itertools.cycle(self.spinner)
        self.color = self._get_color_code(color)
        self.text_color = "\033[37m"
        self.reset_color = "\033[0m"
        self.running = False
        self.thread = None
        self.lock = threading.Lock()
        self.is_paused = False

    def _get_color_code(self, color):
        colors = {
            "black": "\033[30m",
            "red": "\033[31m",
            "green": "\033[32m",
            "yellow": "\033[33m",
            "blue": "\033[34m",
            "magenta": "\033[35m",
            "cyan": "\033[36m",
            "white": "\033[37m",
        }
        return "\033[1m" + colors.get(color.lower(), "\033[37m")

    def _start(self):
        self.running = True
        self.thread = threading.Thread(target=self._animate, daemon=True)
        self.thread.start()

    def _stop(self, success=True):
        with self.lock:
            self.running = False
            if self.thread is not None:
                self.thread.join()
            sys.stdout.write("\r")
            sys.stdout.flush()
            status_icon = emoji.emojize("✅") if success else emoji.emojize("❌")
            sys.stdout.write(
                f"{self.color}{status_icon}{self.reset_color} {self.text_color}{self.text} Done!{self.reset_color}\n"
            )
            sys.stdout.flush()

    def _pause(self):
        with self.lock:
            self.is_paused = True
            sys.stdout.write("\r")
            sys.stdout.flush()

    def _resume(self):
        with self.lock:
            self.is_paused = False

    def _animate(self):
        while self.running:
            if not self.is_paused:
                with self.lock:
                    sys.stdout.write(
                        f"\r{self.color}{next(self.spinner)}{self.reset_color} {self.text_color}{self.text}{self.reset_color} "
                    )
                    sys.stdout.flush()
            time.sleep(0.1)


def show_loader(text="\nLoading...\n", color="cyan"):
    def decorator(func):
        def wrapper(*args, **kwargs):
            loader = Spinner(text=text, color=color)
            loader._start()
            try:
                result = func(*args, **kwargs)
                loader._stop(success=True)
                return result
            except Exception as e:
                loader._stop(success=False)
                raise e

        return wrapper

    return decorator
