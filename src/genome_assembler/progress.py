import shutil
import sys
import time

class ProgressBar:
    def __init__(self, max_value: int) -> None:
        self._current = 0
        self.max_value = max_value

    def update(self, value: int) -> None:
        current = int(value / self.max_value * 100)
        if current != self._current:
            self._current = current
            self.display()

    def display(self) -> None:
        columns, _ = shutil.get_terminal_size()
        sys.stdout.write(f"\r[%-{int(0.8 * columns)}s] %d%%" % ("#" * int(self._current * 0.8 * columns / 100), self._current))
        sys.stdout.flush()

    def finish(self) -> None:
        self.update(self.max_value)
        print()

def main() -> None:
    bar = ProgressBar(50)
    for i in range(50):
        time.sleep(0.1)
        bar.update(i)
    bar.finish()

if __name__ == "__main__":
    main()
