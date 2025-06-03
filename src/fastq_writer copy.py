import threading
import queue
import shutil
import time

# Create a thread-safe queue to hold (temp1, temp2) tuples
file_queue = queue.Queue()

# Set your final output files
final_file_1 = "final_output_1.txt"
final_file_2 = "final_output_2.txt"

# Worker function
def append_worker():
    while True:
        item = file_queue.get()
        if item is None:
            break  # Exit signal received

        temp1, temp2 = item
        try:
            with open(final_file_1, 'ab') as out1, open(temp1, 'rb') as in1:
                shutil.copyfileobj(in1, out1)
            with open(final_file_2, 'ab') as out2, open(temp2, 'rb') as in2:
                shutil.copyfileobj(in2, out2)
        except Exception as e:
            print(f"Error processing {temp1}, {temp2}: {e}")
        finally:
            file_queue.task_done()

# Start the worker thread
worker = threading.Thread(target=append_worker)
worker.start()

# Example main processing loop
for i in range(5):
    # Simulate creation of temp files
    temp1 = f"temp1_{i}.txt"
    temp2 = f"temp2_{i}.txt"
    
    with open(temp1, 'w') as f:
        f.write(f"Data from temp1 file {i}\n")
    with open(temp2, 'w') as f:
        f.write(f"Data from temp2 file {i}\n")

    file_queue.put((temp1, temp2))  # Send file pair to worker

# Signal the worker to shut down and wait for completion
file_queue.put(None)
worker.join()