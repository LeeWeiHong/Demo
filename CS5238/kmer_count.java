import java.io.*;
import java.util.*;

public class kmer_count {
	static final int MEMORY_SIZE = 1000000, HASHSIZE = 250007;
	static final double LOAD_FACTOR = 0.7;
	static final boolean RUNTIME_FLAG = false;
	static final PrintWriter out = new PrintWriter(new BufferedOutputStream(System.out));

	static String filename;
	static int kmerLength, frequency, kmerLimit, numberOfKmers, jellyfishSize;

	// Naive Bloom Filter variables
	static int totalIndices, hashFactor, leftoverKmers, fileNumber = 0;
	static Set<Integer> kmerFrequencySet;

	public static void main(String[] args) throws Exception {
		long startTime = System.currentTimeMillis();

		checkInputParameters(args);
		computeKmerLimit();
		readNumberOfKmers();

		if (numberOfKmers <= kmerLimit) {
			jellyfishSize = numberOfKmers;
			countKmers();
		} else {
			jellyfishSize = (int) Math.ceil(numberOfKmers / kmerLimit);

			countKmersFrequency();

			if (totalIndices < 100) {
				countKmersWithFrequencySet();
			} else {
				for (int i = 1;; i++) {
					hashFactor = i;
					int previousTotalIndices = totalIndices;
					countKmersFrequencyAgain();
					if (totalIndices < 100 || totalIndices >= previousTotalIndices) {
						break;
					}
				}

				countKmersWithWrittenOutput();
			}
		}

		if (RUNTIME_FLAG) {
			long endTime = System.currentTimeMillis();
			out.println("Total runtime = " + (endTime - startTime) + " milliseconds.");
		}

		out.close();
	}

	static void checkInputParameters(String[] args) {
		if (args.length != 3) {
			out.println("Incorrect number of parameters entered.");
			out.println("Correct input format: <filename> <kmer length> <frequency>");
			out.close();
			System.exit(0);
		} else {
			filename = args[0];
			kmerLength = Integer.parseInt(args[1]);
			frequency = Integer.parseInt(args[2]);
		}
	}

	static void computeKmerLimit() {
		kmerLimit = (int) ((LOAD_FACTOR * MEMORY_SIZE) / (2 * kmerLength + 32));
	}

	static void readNumberOfKmers() throws Exception {
		InputStream inputStream = null;

		try {
			inputStream = new FileInputStream(filename);
		} catch (FileNotFoundException e) {
			out.println("Fasta file not found.");
			out.close();
			System.exit(0);
		}

		BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

		char[] kmerArray = new char[kmerLength];
		numberOfKmers = 0;

		while (true) {
			String sequenceName = reader.readLine();
			if (sequenceName == null) {
				break;
			}

			int read, index = 0;
			while ((read = reader.read()) != -1) {
				char c = Character.toLowerCase((char) read);

				if (c == '\n') {
					break;
				}

				if (c == '\r') {
					continue;
				}

				if (c == ' ') {
					index = 0;
					continue;
				}

				if (index != kmerLength) {
					kmerArray[index++] = c;
				} else {
					for (int i = 0; i < kmerLength - 1; i++) {
						kmerArray[i] = kmerArray[i + 1];
					}

					kmerArray[kmerLength - 1] = c;
				}

				if (index == kmerLength && doesNotContainsError(kmerArray)) {
					numberOfKmers++;
				}
			}
		}

		inputStream.close();
		reader.close();
	}

	static void countKmers() throws Exception {
		InputStream inputStream = new FileInputStream(filename);
		BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

		Jellyfish jellyfish = new Jellyfish(jellyfishSize, LOAD_FACTOR);
		char[] kmerArray = new char[kmerLength];
		char[] reverseComplementKmerArray = new char[kmerLength];

		while (true) {
			String sequenceName = reader.readLine();
			if (sequenceName == null) {
				break;
			}

			int read, index = 0;
			while ((read = reader.read()) != -1) {
				char c = Character.toLowerCase((char) read);

				if (c == '\n') {
					break;
				}

				if (c == '\r') {
					continue;
				}

				if (c == ' ') {
					index = 0;
					continue;
				}

				if (index != kmerLength) {
					kmerArray[index++] = c;
				} else {
					for (int i = 0; i < kmerLength - 1; i++) {
						kmerArray[i] = kmerArray[i + 1];
					}

					kmerArray[kmerLength - 1] = c;
				}

				if (index == kmerLength && doesNotContainsError(kmerArray)) {
					if (checkReverseComplement(kmerArray, reverseComplementKmerArray)) {
						jellyfish.addKmer(kmerArray);
					} else {
						jellyfish.addKmer(reverseComplementKmerArray);
					}
				}
			}
		}

		inputStream.close();
		reader.close();
		jellyfish.printKmerCount(frequency, out);
	}

	static void countKmersFrequency() throws Exception {
		InputStream inputStream = new FileInputStream(filename);
		BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

		short[] kmerFrequency = new short[HASHSIZE];
		char[] kmerArray = new char[kmerLength];
		char[] reverseComplementKmerArray = new char[kmerLength];

		while (true) {
			String sequenceName = reader.readLine();
			if (sequenceName == null) {
				break;
			}

			int read, index = 0;
			while ((read = reader.read()) != -1) {
				char c = Character.toLowerCase((char) read);

				if (c == '\n') {
					break;
				}

				if (c == '\r') {
					continue;
				}

				if (c == ' ') {
					index = 0;
					continue;
				}

				if (index != kmerLength) {
					kmerArray[index++] = c;
				} else {
					for (int i = 0; i < kmerLength - 1; i++) {
						kmerArray[i] = kmerArray[i + 1];
					}

					kmerArray[kmerLength - 1] = c;
				}

				if (index == kmerLength && doesNotContainsError(kmerArray)) {
					int hashsum = 0;

					if (checkReverseComplement(kmerArray, reverseComplementKmerArray)) {
						hashsum = hashForKmerFrequency(kmerArray);
					} else {
						hashsum = hashForKmerFrequency(reverseComplementKmerArray);
					}

					kmerFrequency[hashsum]++;
				}
			}
		}

		inputStream.close();
		reader.close();

		totalIndices = 0;
		for (int i = 0; i < kmerFrequency.length; i++) {
			if (kmerFrequency[i] >= frequency) {
				totalIndices++;
			}
		}
		
		if (totalIndices < 100) {
			BufferedWriter writer = new BufferedWriter(new FileWriter("kmersFrequencyIndex.txt"));
			for (int i = 0; i < kmerFrequency.length; i++) {
				if (kmerFrequency[i] >= frequency) {
					writer.write("" + i);
					writer.newLine();
				}
			}

			writer.close();
			readKmersFrequency(totalIndices);
		} else {
			writeKmersWithFrequency(kmerFrequency);
		}
	}

	static void readKmersFrequency(int totalKmers) throws Exception {
		InputStream inputStream = new FileInputStream("kmersFrequencyIndex.txt");
		BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

		kmerFrequencySet = new HashSet<Integer>();

		while (totalKmers-- > 0) {
			kmerFrequencySet.add(Integer.parseInt(reader.readLine()));
		}

		inputStream.close();
		reader.close();
	}

	static void countKmersWithFrequencySet() throws Exception {
		InputStream inputStream = new FileInputStream(filename);
		BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
		BufferedWriter writer = new BufferedWriter(new FileWriter("kmers" + fileNumber + ".txt"));

		Jellyfish jellyfish = new Jellyfish(jellyfishSize, LOAD_FACTOR);
		char[] kmerArray = new char[kmerLength];
		char[] reverseComplementKmerArray = new char[kmerLength];
		leftoverKmers = 0;
		boolean useKmerArray = true, kmerLimitReached = false;
		int hashsum = 0;

		while (true) {
			String sequenceName = reader.readLine();
			if (sequenceName == null) {
				break;
			}

			int read, index = 0;
			while ((read = reader.read()) != -1) {
				char c = Character.toLowerCase((char) read);

				if (c == '\n') {
					break;
				}

				if (c == '\r') {
					continue;
				}

				if (c == ' ') {
					index = 0;
					continue;
				}

				if (index != kmerLength) {
					kmerArray[index++] = c;
				} else {
					for (int i = 0; i < kmerLength - 1; i++) {
						kmerArray[i] = kmerArray[i + 1];
					}

					kmerArray[kmerLength - 1] = c;
				}

				if (index == kmerLength && doesNotContainsError(kmerArray)) {
					useKmerArray = checkReverseComplement(kmerArray, reverseComplementKmerArray);

					if (useKmerArray) {
						hashsum = hashForKmerFrequency(kmerArray);
					} else {
						hashsum = hashForKmerFrequency(reverseComplementKmerArray);
					}

					if (kmerFrequencySet.contains(hashsum)) {
						if (useKmerArray) {
							kmerLimitReached = jellyfish.addKmer(kmerArray);
						} else {
							kmerLimitReached = jellyfish.addKmer(reverseComplementKmerArray);
						}

						if (kmerLimitReached) {
							leftoverKmers++;
							if (useKmerArray) {
								writer.write(kmerArray, 0, kmerLength);
							} else {
								writer.write(reverseComplementKmerArray, 0, kmerLength);
							}
							writer.newLine();
						}
					}
				}
			}
		}

		inputStream.close();
		reader.close();
		writer.close();
		
		jellyfish.printKmerCount(frequency, out);

		if (kmerLimitReached) {
			countKmersWithWrittenOutput();
		}
	}

	static void writeKmersWithFrequency(short[] kmerFrequency) throws Exception {
		InputStream inputStream = new FileInputStream(filename);
		BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
		BufferedWriter writer = new BufferedWriter(new FileWriter("kmers" + fileNumber + ".txt"));

		char[] kmerArray = new char[kmerLength];
		char[] reverseComplementKmerArray = new char[kmerLength];
		leftoverKmers = 0;

		while (true) {
			String sequenceName = reader.readLine();
			if (sequenceName == null) {
				break;
			}

			int read, index = 0;
			while ((read = reader.read()) != -1) {
				char c = Character.toLowerCase((char) read);

				if (c == '\n') {
					break;
				}

				if (c == '\r') {
					continue;
				}

				if (c == ' ') {
					index = 0;
					continue;
				}

				if (index != kmerLength) {
					kmerArray[index++] = c;
				} else {
					for (int i = 0; i < kmerLength - 1; i++) {
						kmerArray[i] = kmerArray[i + 1];
					}

					kmerArray[kmerLength - 1] = c;
				}

				if (index == kmerLength && doesNotContainsError(kmerArray)) {
					int hashsum = 0;
					boolean useKmerArray = checkReverseComplement(kmerArray, reverseComplementKmerArray);

					if (useKmerArray) {
						hashsum = hashForKmerFrequency(kmerArray);
					} else {
						hashsum = hashForKmerFrequency(reverseComplementKmerArray);
					}

					if (kmerFrequency[hashsum] >= frequency) {
						leftoverKmers++;
						if (useKmerArray) {
							writer.write(kmerArray, 0, kmerLength);
						} else {
							writer.write(reverseComplementKmerArray, 0, kmerLength);
						}
						writer.newLine();
					}
				}
			}
		}

		inputStream.close();
		reader.close();
		writer.close();
	}

	static void countKmersFrequencyAgain() throws Exception {
		InputStream inputStream = new FileInputStream("kmers" + fileNumber + ".txt");
		BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

		short[] kmerFrequency = new short[HASHSIZE];
		char[] kmerArray = new char[kmerLength];
		char[] reverseComplementKmerArray = new char[kmerLength];
		int numberOfReads = leftoverKmers, hashsum = 0;

		while (numberOfReads-- > 0) {
			reader.read(kmerArray, 0, kmerLength);
			reader.read();
			reader.read();

			if (checkReverseComplement(kmerArray, reverseComplementKmerArray)) {
				hashsum = hashForKmerFrequency(kmerArray);
			} else {
				hashsum = hashForKmerFrequency(reverseComplementKmerArray);
			}

			kmerFrequency[hashsum]++;
		}

		inputStream.close();
		reader.close();

		totalIndices = 0;
		for (int i = 0; i < kmerFrequency.length; i++) {
			if (kmerFrequency[i] >= frequency) {
				totalIndices++;
			}
		}

		if (totalIndices >= 100) {
			writeKmersWithFrequencyAgain(kmerFrequency);
		}
	}

	static void writeKmersWithFrequencyAgain(short[] kmerFrequency) throws Exception {
		InputStream inputStream = new FileInputStream("kmers" + fileNumber++ + ".txt");
		BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
		BufferedWriter writer = new BufferedWriter(new FileWriter("kmers" + fileNumber + ".txt"));

		char[] kmerArray = new char[kmerLength];
		char[] reverseComplementKmerArray = new char[kmerLength];
		boolean useKmerArray = true;
		int numberOfReads = leftoverKmers, hashsum = 0;
		leftoverKmers = 0;

		while (numberOfReads-- > 0) {
			reader.read(kmerArray, 0, kmerLength);
			reader.read();
			reader.read();

			useKmerArray = checkReverseComplement(kmerArray, reverseComplementKmerArray);

			if (useKmerArray) {
				hashsum = hashForKmerFrequency(kmerArray);
			} else {
				hashsum = hashForKmerFrequency(reverseComplementKmerArray);
			}

			if (kmerFrequency[hashsum] >= frequency) {
				leftoverKmers++;
				if (useKmerArray) {
					writer.write(kmerArray, 0, kmerLength);
				} else {
					writer.write(reverseComplementKmerArray, 0, kmerLength);
				}
				writer.newLine();
			}
		}

		inputStream.close();
		reader.close();
		writer.close();
	}

	static void countKmersWithWrittenOutput() throws Exception {
		while (leftoverKmers != 0) {
			InputStream inputStream = new FileInputStream("kmers" + fileNumber++ + ".txt");
			BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
			BufferedWriter writer = new BufferedWriter(new FileWriter("kmers" + fileNumber + ".txt"));

			int jellyfishSize = (int) Math.ceil(numberOfKmers / kmerLimit);
			Jellyfish jellyfish = new Jellyfish(jellyfishSize, LOAD_FACTOR);
			char[] kmerArray = new char[kmerLength];
			char[] reverseComplementKmerArray = new char[kmerLength];
			boolean useKmerArray, kmerLimitReached = false;
			int numberOfReads = leftoverKmers;
			leftoverKmers = 0;

			while (numberOfReads-- > 0) {
				reader.read(kmerArray, 0, kmerLength);
				reader.read();
				reader.read();
				useKmerArray = checkReverseComplement(kmerArray, reverseComplementKmerArray);

				if (useKmerArray) {
					kmerLimitReached = jellyfish.addKmer(kmerArray);
				} else {
					kmerLimitReached = jellyfish.addKmer(reverseComplementKmerArray);
				}

				if (kmerLimitReached) {
					leftoverKmers++;
					if (useKmerArray) {
						writer.write(kmerArray, 0, kmerLength);
					} else {
						writer.write(reverseComplementKmerArray, 0, kmerLength);
					}
					writer.newLine();
				}
			}

			inputStream.close();
			reader.close();
			writer.close();

			jellyfish.printKmerCount(frequency, out);
		}
	}

	static boolean doesNotContainsError(char[] kmerArray) {
		for (char c : kmerArray) {
			if (c == 'n') {
				return false;
			}
		}

		return true;
	}

	static boolean checkReverseComplement(char[] kmerArray, char[] reverseComplementKmerArray) {
		for (int i = 0, j = kmerArray.length - 1; i < kmerArray.length; i++, j--) {
			if (kmerArray[i] == 'a') {
				reverseComplementKmerArray[j] = 't';
			} else if (kmerArray[i] == 't') {
				reverseComplementKmerArray[j] = 'a';
			} else if (kmerArray[i] == 'c') {
				reverseComplementKmerArray[j] = 'g';
			} else if (kmerArray[i] == 'g') {
				reverseComplementKmerArray[j] = 'c';
			}
		}

		for (int i = 0; i < kmerArray.length; i++) {
			if (kmerArray[i] < reverseComplementKmerArray[i]) {
				break;
			} else if (kmerArray[i] > reverseComplementKmerArray[i]) {
				return false;
			}
		}

		return true;
	}

	static int hashForKmerFrequency(char[] kmerArray) {
		int hashsum = 0, base = 0, firstBase = 0;

		if (kmerArray[0] == 'a') {
			hashsum = 0 + hashFactor;
		} else if (kmerArray[0] == 'c') {
			hashsum = 1 + hashFactor;
		} else if (kmerArray[0] == 'g') {
			hashsum = 2 + hashFactor;
		} else if (kmerArray[0] == 't') {
			hashsum = 3 + hashFactor;
		}

		firstBase = hashsum %= 4;

		for (int i = 1; i < kmerArray.length; i++) {
			if (kmerArray[i] == 'a') {
				base = 0 + hashFactor;
			} else if (kmerArray[i] == 'c') {
				base = 1 + hashFactor;
			} else if (kmerArray[i] == 'g') {
				base = 2 + hashFactor;
			} else if (kmerArray[i] == 't') {
				base = 3 + hashFactor;
			}

			base %= 4;
			hashsum = (hashsum * 4 + base) % HASHSIZE;
		}

		return hashsum = (hashsum >> firstBase) % HASHSIZE;
	}
}

class Jellyfish {
	private String[] hashTable;
	private int[] countTable;
	private int hashsize, numberOfKmers, kmerLimit;

	public Jellyfish(int hashsize, double loadFactor) {
		this.hashsize = (int) (hashsize / loadFactor);
		hashTable = new String[this.hashsize];
		countTable = new int[this.hashsize];
		numberOfKmers = 0;
		kmerLimit = hashsize;
	}

	boolean addKmer(char[] kmerArray) {
		int hashsum = hash(kmerArray);

		if (hashsum != -1) {
			hashTable[hashsum] = new String(kmerArray);
			countTable[hashsum]++;
			return false;
		}

		return true;
	}

	void printKmerCount(int frequency, PrintWriter out) {
		for (int i = 0; i < hashsize; i++) {
			if (hashTable[i] != null && countTable[i] >= frequency) {
				out.println(countTable[i] + " " + hashTable[i]);
			}
		}
	}

	int hash(char[] kmerArray) {
		int hashsum = 0, base = 0;

		if (kmerArray[0] == 'c') {
			hashsum = 1;
		} else if (kmerArray[0] == 'g') {
			hashsum = 2;
		} else if (kmerArray[0] == 't') {
			hashsum = 3;
		}

		for (int i = 1; i < kmerArray.length; i++) {
			if (kmerArray[i] == 'a') {
				base = 0;
			} else if (kmerArray[i] == 'c') {
				base = 1;
			} else if (kmerArray[i] == 'g') {
				base = 2;
			} else if (kmerArray[i] == 't') {
				base = 3;
			}

			hashsum = (hashsum * 4 + base) % hashsize;
		}

		hashsum = hashsum % hashsize;

		String kmer = new String(kmerArray);
		while (hashTable[hashsum] != null && !hashTable[hashsum].equals(kmer)) {
			hashsum = (hashsum + 1) % hashsize;
		}

		if (numberOfKmers >= kmerLimit && hashTable[hashsum] == null) {
			return -1;
		}

		if (hashTable[hashsum] == null) {
			numberOfKmers++;
		}

		return hashsum;
	}
}