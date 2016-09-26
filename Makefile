all:
	gcc -std=gnu99 -pedantic -lm -o build main.c p2random.c tree.c -lrt -O2
