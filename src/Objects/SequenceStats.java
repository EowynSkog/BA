package Objects;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import inb_bach.MainClass;

public class SequenceStats {
	private String Sequence;
	// Array Params: Base1. Base2, Base3, x, 0/1
	// x=Base in front or after
	// 0= Base in Front of Triplet
	// 1=Base after triplet
	private int[][][][][] rawData = new int[5][5][5][5][2];
	private int[][][][][][][] rawDataTriplet = new int[5][5][5][5][5][5][2];
	private String[] Bases = { "T", "C", "A", "G", "N" };
	static DecimalFormat df = new DecimalFormat("0.0000"); // ?
	private double[] Base_aPriori;
	private int[] BaseCount;
	private double[][][] Triplet_aPriori;
	private int[][][] TripletCount;
	private double[][][][][] TripletTransition;
	private double[][][][][][][] TripletCompleteTransition;
	private double[][] BaseTransition;
	private int[][] BaseTransitionCount;

	public SequenceStats(String Sequence) {
		this.Sequence = Sequence;
		System.out.println("Calculatinng raw Statistics-Matrix..");
		// calculateRawDataTripletTransiton();
		// calculateTripletCompleteTransition();
		calculateRawData(); // Zählt Triplettvorkommen mit base davor und base
		// // danach
		 System.out.println("Processing raw Data...");
		 calculateBase_aPriori();
		 calculateTriplet_aPriori();
		 calculateBaseTransition();
		calculateTripletTransition();
	}

	// Calculates NA weightings (Nucleotide a priori probability) /Zählt die in
	// BaseCount an Stelle 1: T, 2: C, 3: A, 4: G, 5:ges Basen Anzahl
	private void calculateBase_aPriori() {
		BaseCount = new int[6];
		int overallCount = 0;
		Base_aPriori = new double[5];
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					for (int l = 0; l < 5; l++) {
						// Zählt einzelne Basen der Tripletts (unabhängig vom
						// Nachfolger
						// Auch die letzen 3 Basen! Denn default Nachfolger als
						// N behandelt)
						BaseCount[i] += rawData[i][j][k][l][1];
						BaseCount[j] += rawData[i][j][k][l][1];
						BaseCount[k] += rawData[i][j][k][l][1];
						overallCount += rawData[i][j][k][l][1] * 3;
					}
				}
			}
		}
		BaseCount[5] = overallCount;
		for (int i = 0; i < 5; i++) {
			Base_aPriori[i] = ((double) BaseCount[i] / (double) overallCount) * 4; // ohne N -> *4
			// warum mal 4?

		}

		System.out.println(
				"Base_aPriori Average: " + (Base_aPriori[0] + Base_aPriori[1] + Base_aPriori[2] + Base_aPriori[3]) / 4); // N-> /5, ohne-> /4
	}

	// Calculates TA weightings
	private void calculateTriplet_aPriori() {
		TripletCount = new int[5][5][5];
		int overallCount = 0;
		Triplet_aPriori = new double[5][5][5];
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					for (int l = 0; l < 5; l++) {

						// String Pattern=Bases[i]+Bases[j]+Bases[k];
						// if
						// (Pattern.equalsIgnoreCase("TAA")||Pattern.equalsIgnoreCase("TAG")||Pattern.equalsIgnoreCase("TGA"))continue;
						// // Stop codons werden nicht mitgezählt
						TripletCount[i][j][k] += rawData[i][j][k][l][1];
						overallCount += rawData[i][j][k][l][1];
					}
				}
			}
		}
		double sum = 0;
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					Triplet_aPriori[i][j][k] = ((double) TripletCount[i][j][k] * 64) / (double) overallCount; // statt 61  -> 64  weil die  Stopcodons mitgezählt  werden
					//wegen N Behandlung *125 statt *64
					sum += Triplet_aPriori[i][j][k];
				}
			}
		}
		sum = sum / 64; //N-> 125, ohne -> 64
		System.out.println("Triplet_aPriori Sum: " + sum);

	}
	// Seves stuff I don't want to crash
	// private void calculateTriplet_aPriori(){
	// int[][][] TripletCount=new int[4][4][4];
	// int overallCount=0;
	// Triplet_aPriori=new double[4][4][4];
	// for(int i=0;i<4;i++){
	// for(int j=0;j<4;j++){
	// for(int k=0;k<4;k++){
	// for(int l=0;l<4;l++){
	// //This will ignore 4 Bases at the beginning - shit happens :P
	// String Pattern=Bases[i]+Bases[j]+Bases[k];
	// if
	// (Pattern.equalsIgnoreCase("TAA")||Pattern.equalsIgnoreCase("TAG")||Pattern.equalsIgnoreCase("TGA"))continue;
	// // Stop codons werden nicht mitgezählt
	// TripletCount[i][j][k]+=rawData[i][j][k][l][1];// Achtung alle Möglichen
	// Kombinationen werden hier als Triplett gezählt!
	// overallCount+=rawData[i][j][k][l][1]; //Anzahl der Tripletts ist mehr als
	// Tripletts Anzahl der Basen/3 ; Stopcodons werden hier nicht mitgezählt
	// }
	// }
	// }
	// }
	// double sum=0;
	// for(int i=0;i<4;i++){
	// for(int j=0;j<4;j++){
	// for(int k=0;k<4;k++){
	// Triplet_aPriori[i][j][k]=((double)TripletCount[i][j][k]*61)/(double)overallCount;
	// //(Anzahl eines Triplets*Anzahle aller Möglichen Unterschiedlichen
	// (=61))/Anzahl aller
	// sum+=Triplet_aPriori[i][j][k];
	// }
	// }
	// }
	// sum=sum/61;
	// System.out.println("Triplet_aPriori Average: "+sum);
	// }

	// activate Seves stuff hear

	// Calculates NT weightings
	private void calculateBaseTransition() {
		BaseTransitionCount = new int[5][5];
		int overallCount = 0;
		BaseTransition = new double[5][5];
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					for (int l = 0; l < 5; l++) {
						BaseTransitionCount[i][j] += rawData[i][j][k][l][1];
						BaseTransitionCount[j][k] += rawData[i][j][k][l][1];
						BaseTransitionCount[k][l] += rawData[i][j][k][l][1];
						overallCount += rawData[i][j][k][l][1] * 3;
					}
				}
			}
		}
		double sum = 0;
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				BaseTransition[i][j] = ((double) BaseTransitionCount[i][j] / (double) overallCount) * 16; // mit 25, ohne N -> *16
				sum += BaseTransition[i][j];
			}
		}
		sum = sum / 16; //mit 25, ohne N-> /16
		System.out.println("BaseTransition Average: " + sum);
	}

	// Calculates TT weightings
	private void calculateTripletTransition() {
		int overallCountFront = 0;
		int overallCountAfter = 0;
		TripletTransition = new double[5][5][5][5][2];
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					for (int l = 0; l < 5; l++) {
						overallCountFront += rawData[i][j][k][l][0];
						overallCountAfter += rawData[i][j][k][l][1];
					}
				}
			}
		}
		double sumFront = 0;
		double sumAfter = 0;
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					for (int l = 0; l < 5; l++) {
						TripletTransition[i][j][k][l][0] = ((double) rawData[i][j][k][l][0]
								/ (double) overallCountFront) * 256;
						TripletTransition[i][j][k][l][1] = ((double) rawData[i][j][k][l][1]
								/ (double) overallCountAfter) * 256; //ohne N-> 256=4^3*4; mit N 625=5^3*5
						sumFront += TripletTransition[i][j][k][l][0];
						sumAfter += TripletTransition[i][j][k][l][1];
					}
				}
			}
		}
		sumFront = sumFront / 256; // ohne N -> 256
		sumAfter = sumAfter / 256; // ohne N -> 256
		System.out.println("Triplet Transition Average: Front:" + sumFront + " After:" + sumAfter);
	}

	private void calculateRawData() {
		// Iterate over whole Sequence;
		for (int i = 0; i < Sequence.length() - 2; i = i + 3) { // im 3er
																// Schritt
																// 1-letzte Base
			int a = -1;
			int b = -1;
			int c = -1;
			int front = -1;
			int after = -1;

			// First Base of Triplet
			switch (Sequence.charAt(i)) {
			case 'T':
				a = 0;
				break;
			case 'C':
				a = 1;
				break;
			case 'A':
				a = 2;
				break;
			case 'G':
				a = 3;
				break;
			case 'N':
				a = 4;
				break;
			default:
				System.out.println("Unerwarteter Buchstabe in der Sequenz");
				continue;
			}

			// Second Base of Triplet
			switch (Sequence.charAt(i + 1)) {
			case 'T':
				b = 0;
				break;
			case 'C':
				b = 1;
				break;
			case 'A':
				b = 2;
				break;
			case 'G':
				b = 3;
				break;
			case 'N':
				b = 4;
				break;
			default:
				System.out.println("Unerwarteter Buchstabe in der Sequenz");
				continue;
			}

			// Third Base of Triplet
			switch (Sequence.charAt(i + 2)) {
			case 'T':
				c = 0;
				break;
			case 'C':
				c = 1;
				break;
			case 'A':
				c = 2;
				break;
			case 'G':
				c = 3;
				break;
			case 'N':
				c = 4;
				break;
			default:
				System.out.println("Unerwarteter Buchstabe in der Sequenz");
				continue;
			}

			// Base After Triplet
			if ((i + 3) == (Sequence.length())) { // speichert default mäßig die
													// Nachfolgebase nach der
													// Sequenz als N ab
				System.out.println("Sequenzende erreicht. Letze Base " + Sequence.charAt(i + 2));
				after = 4;
			} else if ((i + 3) < Sequence.length()) {
				switch (Sequence.charAt(i + 3)) {
				case 'T':
					after = 0;
					break;
				case 'C':
					after = 1;
					break;
				case 'A':
					after = 2;
					break;
				case 'G':
					after = 3;
					break;
				case 'N':
					after = 4;
					break;
				}

			}

			// Base Before Triplet
			if ((i - 1) >= 0) {
				switch (Sequence.charAt(i - 1)) {
				case 'T':
					front = 0;
					break;
				case 'C':
					front = 1;
					break;
				case 'A':
					front = 2;
					break;
				case 'G':
					front = 3;
					break;
				case 'N':
					front = 4;
					break;
				}
			}
			// Ignoring Special Characters
			if (front == -1 && after == -1)
				continue;

			if (front != -1) {
				rawData[a][b][c][front][0]++;
			}
			if (after != -1) {
				rawData[a][b][c][after][1]++;
			}
		}
	}

	private void calculateRawDataTripletTransiton() {
		// Iterate over whole Sequence;
		for (int i = 0; i < Sequence.length() - 2; i = i + 3) { // im 3er
																// Schritt

			int a = -1;
			int b = -1;
			int c = -1;
			int front = -1;
			int fb = -1;
			int fc = -1;
			int after = -1;
			int ab = -1;
			int ac = -1;

			// First Base of Triplet
			switch (Sequence.charAt(i)) {
			case 'T':
				a = 0;
				break;
			case 'C':
				a = 1;
				break;
			case 'A':
				a = 2;
				break;
			case 'G':
				a = 3;
				break;
			case 'N':
				a = 4;
				break;
			default:
				continue;
			}

			// Second Base of Triplet
			switch (Sequence.charAt(i + 1)) {
			case 'T':
				b = 0;
				break;
			case 'C':
				b = 1;
				break;
			case 'A':
				b = 2;
				break;
			case 'G':
				b = 3;
				break;
			case 'N':
				b = 4;
				break;
			default:
				continue;
			}

			// Third Base of Triplet
			switch (Sequence.charAt(i + 2)) {
			case 'T':
				c = 0;
				break;
			case 'C':
				c = 1;
				break;
			case 'A':
				c = 2;
				break;
			case 'G':
				c = 3;
				break;
			case 'N':
				c = 4;
				break;
			default:
				continue;
			}

			// Base After Triplet
			if ((i + 3) < Sequence.length()) {
				switch (Sequence.charAt(i + 3)) {
				case 'T':
					after = 0;
					break;
				case 'C':
					after = 1;
					break;
				case 'A':
					after = 2;
					break;
				case 'G':
					after = 3;
					break;
				case 'N':
					after = 4;
					break;
				}

			}
			if ((i + 4) < Sequence.length()) {
				switch (Sequence.charAt(i + 4)) {
				case 'T':
					ab = 0;
					break;
				case 'C':
					ab = 1;
					break;
				case 'A':
					ab = 2;
					break;
				case 'G':
					ab = 3;
					break;
				case 'N':
					ab = 4;
					break;
				}

			}
			if ((i + 5) < Sequence.length()) {
				switch (Sequence.charAt(i + 5)) {
				case 'T':
					ac = 0;
					break;
				case 'C':
					ac = 1;
					break;
				case 'A':
					ac = 2;
					break;
				case 'G':
					ac = 3;
					break;
				case 'N':
					ac = 4;
					break;
				}

			}

			// Base Before Triplet
			if ((i - 1) >= 0) {
				switch (Sequence.charAt(i - 1)) {
				case 'T':
					front = 0;
					break;
				case 'C':
					front = 1;
					break;
				case 'A':
					front = 2;
					break;
				case 'G':
					front = 3;
					break;
				case 'N':
					front = 4;
					break;
				}
			}
			if ((i - 2) >= 0) {
				switch (Sequence.charAt(i - 2)) {
				case 'T':
					fb = 0;
					break;
				case 'C':
					fb = 1;
					break;
				case 'A':
					fb = 2;
					break;
				case 'G':
					fb = 3;
					break;
				case 'N':
					fb = 4;
					break;
				}
			}
			if ((i - 3) >= 0) {
				switch (Sequence.charAt(i - 3)) {
				case 'T':
					fc = 0;
					break;
				case 'C':
					fc = 1;
					break;
				case 'A':
					fc = 2;
					break;
				case 'G':
					fc = 3;
					break;
				case 'N':
					fc = 4;
					break;
				}
			}
			// Ignoring Special Characters
			if (front == -1 && ac == -1)
				continue;

			if (fc != -1) {
				rawDataTriplet[a][b][c][fc][fb][front][0]++;
			}
			if (ac != -1) {
				rawDataTriplet[a][b][c][after][ab][ac][1]++;
			}
		}
	}

	private void calculateTripletCompleteTransition() {
		int overallCountFront = 0;
		int overallCountAfter = 0;
		TripletCompleteTransition = new double[5][5][5][5][5][5][2];
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				for (int k = 0; k < 5; k++) {
					for (int l = 0; l < 5; l++) {
						for (int m = 0; m < 5; m++) {
							for (int n = 0; n < 5; n++) {
								TripletCompleteTransition[i][j][k][l][m][n][0] = rawDataTriplet[i][j][k][l][m][n][0];
								TripletCompleteTransition[i][j][k][l][m][n][1] = rawDataTriplet[i][j][k][l][m][n][1];
								overallCountFront += rawDataTriplet[i][j][k][l][m][n][0];
								overallCountAfter += rawDataTriplet[i][j][k][l][m][n][1];
							}
						}
					}
				}
			}
		}
		double sumFront = 0;
		double sumAfter = 0;

	}

	// Difference by Element between matrices //Wozu das denn?
	public double[][] MatrixDiff(double[][] M1, double[][] M2) {
		double[][] Erg = new double[4][4];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Erg[i][j] = M1[i][j] - M2[i][j];
			}
		}
		double sum = 0;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				sum += Math.abs(Erg[i][j]);
			}
		}
		System.out.println("Gesamtdiferenz: " + df.format(sum));
		PrintMatrix(Erg);
		return Erg;
	}

	// Prints the 4x4 transition matrix in the console //???
	public void PrintMatrix(double[][] Proz) {
		System.out.println("Vertikal: s(n) horizontal: s(n+1)");
		System.out.println("--- T ------- C ------- A ------- G");
		System.out.println("T " + df.format(Proz[0][0] / 4) + " -- " + df.format(Proz[1][0] / 4) + " -- "
				+ df.format(Proz[2][0] / 4) + " -- " + df.format(Proz[3][0] / 4));
		System.out.println("C " + df.format(Proz[0][1] / 4) + " -- " + df.format(Proz[1][1] / 4) + " -- "
				+ df.format(Proz[2][1] / 4) + " -- " + df.format(Proz[3][1] / 4));
		System.out.println("A " + df.format(Proz[0][2] / 4) + " -- " + df.format(Proz[1][2] / 4) + " -- "
				+ df.format(Proz[2][2] / 4) + " -- " + df.format(Proz[3][2] / 4));
		System.out.println("G " + df.format(Proz[0][3] / 4) + " -- " + df.format(Proz[1][3] / 4) + " -- "
				+ df.format(Proz[2][3] / 4) + " -- " + df.format(Proz[3][3] / 4));
	}

	public int[][][][][] getRawData() {
		return rawData;
	}

	public double[] getBase_aPriori() {
		return Base_aPriori;
	}

	public int[] getBaseCount() {
		return BaseCount;
	}

	public double[][][] getTriplet_aPriori() {
		return Triplet_aPriori;
	}

	public int[][][] getTripletCount() {
		return TripletCount;
	}

	public double[][][][][] getTripletTransition() {
		return TripletTransition;
	}

	public double[][] getBaseTransition() {
		return BaseTransition;
	}

	public int[][] getBaseTransitionCount() {
		return BaseTransitionCount;
	}

	public double[][][][][][][] getTripletCompleteTransition() {
		return TripletCompleteTransition;
	}

}
