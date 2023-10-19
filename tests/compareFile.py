import sys
import hashlib


def hashfile(file):

	# A arbitrary (but fixed) buffer 
	# size (change accordingly)
	# 65536 = 65536 bytes = 64 kilobytes 
	BUF_SIZE = 65536

	# Initializing the sha256() method
	sha256 = hashlib.sha256()

	# Opening the file provided as
	# the first commandline argument
	with open(file, 'rb') as f:
		
		while True:
			
			# reading data = BUF_SIZE from
			# the file and saving it in a
			# variable
			data = f.read(BUF_SIZE)

			# True if eof = 1
			if not data:
				break
	
			# Passing that data to that sh256 hash
			# function (updating the function with
			# that data)
			sha256.update(data)

	
	# sha256.hexdigest() hashes all the input
	# data passed to the sha256() via sha256.update()
	# Acts as a finalize method, after which
	# all the input data gets hashed hexdigest()
	# hashes the data, and returns the output
	# in hexadecimal format
	return sha256.hexdigest()

# Calling hashfile() function to obtain hashes
# of the files, and saving the result
# in a variable
f1_hash = hashfile(sys.argv[1])
f2_hash = hashfile(sys.argv[2])

# Doing primitive string comparison to 
# check whether the two hashes match or not
if f1_hash == f2_hash:
	print("Both files are same")
	print(f"Hash: {f1_hash}")

else:
	print("Files are different!")
	print(f"Hash of File 1: {f1_hash}")
	print(f"Hash of File 2: {f2_hash}")
