import time

with open("output.txt", "w") as f:
    f.write("Numpy not imported yet.")
try:
    import numpy
    # Example Python code to test
    def main():
        print("Starting the job...")
        with open("output.txt", "w") as f:
            f.write("This is a test file.")
        #     f.write(help(numpy))
        # help(numpy)
        time.sleep(1)
        
        print("Job completed!")
except ImportError as e:
    with open("output.txt", "w") as f:
        f.write(f"Error: {e}")
if __name__ == "__main__":
    main()
