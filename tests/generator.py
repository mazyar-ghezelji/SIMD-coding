import random
seed = ""
bases = ["A", "C", "T", "G"]
random.seed(10)
for i in range(200): # length of each read would be 200 bp 
    seed += (bases[random.randint(0, 3)])

f = open("sequences.txt", "w")
t = open("sequences.fa", "w")
for i in range(10001):
    f.write(">read"+str(i)+"\n")
    t.write(">read"+str(i)+"\n")
    seq = list(seed)
    freq = random.randint(160, 180)
    for j in range(freq):
        seq[random.randint(0, 199)] = bases[random.randint(0, 3)]
    a = "".join(seq)
    f.write(a+"\n")
    t.write(a+"\n")