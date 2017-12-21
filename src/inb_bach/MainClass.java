package inb_bach;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;

import javax.sound.midi.Sequence;

import org.apache.commons.io.FileUtils;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;

import Objects.GeneCode;
import Objects.SequenceStats;

public class MainClass {
	// von Berit eingefügt Runtime Test
	private static final long MEGABYTE = 1024L * 1024L;

	public static long bytesToMegabytes(long bytes) {
		return bytes / MEGABYTE;
	}
	// runtime Test ende

	static String[] Code = { "Leu", "Pro", "His", "Gln", "Arg", "Ile", "Met", "Thr", "Asn", "Lys", "Ser", "Val", "Ala",
			"Asp", "Glu", "Gly", "Phe", "Tyr", "Cys", "Trp" };
	static DecimalFormat df = new DecimalFormat("0.0000");
	static GenBankConnection conn = new GenBankConnection();
	public static double[] baseAprioriWeights;
	public static int[] baseAprioriCount;
	public static double[][][] tripletAprioriWeights;
	public static int[][][] tripletAprioriCount;
	public static double[][] baseTransitionWeights;
	public static int[][] baseTransitionCount;
	public static double[][][][][] tripletTransitionWeights;
	public static int[][][][][] tripletTransitionCount;
	public static double[][][][][][][] tripletComleteTransition;
	public static boolean baseAprioriEnabled = false;
	public static boolean tripletAprioriEnabled = false;
	public static boolean baseTransitionEnabled = false;
	public static boolean tripletTransitionEnabled = false;
	public static int TransitionTransversionBias = 1;

	public static final String sequenceFileName = "human_CDS_biomart";

	public static void main(String[] args) {
		// ###################################################################
		// ######################### DEBUG AREA ##########################
		// ###################################################################
		// statistic Berits try

		// mixed Start
		List<DNASequence> Mixed = conn.LoadMixedFile(sequenceFileName);
		StringBuilder builder = new StringBuilder();
		for (DNASequence seq : Mixed) {
			builder.append(seq.getSequenceAsString());
		}
		String MixedSeq = builder.toString();
		System.out.println("Mixed Length: " + MixedSeq.length());

		SequenceStats Stat = new SequenceStats(MixedSeq);
		// Mixed end
		// Single Sequence start
		// DNASequence MySeq = conn.LoadFastaFile(sequenceFileName);
		// SequenceStats Stat = new SequenceStats(MySeq.getSequenceAsString());

		// statistics
		baseAprioriWeights = Stat.getBase_aPriori();
		// baseAprioriCount = Stat.getBaseCount();
		tripletAprioriWeights = Stat.getTriplet_aPriori();
		// tripletAprioriCount = Stat.getTripletCount();
		baseTransitionWeights = Stat.getBaseTransition();
		// baseTransitionCount = Stat.getBaseTransitionCount();
		tripletTransitionWeights = Stat.getTripletTransition();
		// tripletTransitionCount = Stat.getRawData();
		// tripletComleteTransition = Stat.getTripletCompleteTransition();

		TransitionTransversionBias = 2;
		setWeightings(true, true, true, true);
		GeneCode G = new GeneCode(Code);
		StabilityCalculator S = new StabilityCalculator(G);

		double MS1 = S.get_BaseDeviation(1);// MS1
		double MS2 = S.get_BaseDeviation(2);// MS2
		double MS3 = S.get_BaseDeviation(3);// MS3
		double MS0 = S.getMS0(MS1, MS2, MS3);// MS0
		double rMS = S.get_ShiftDeviation(1);// rMS
		double lMS = S.get_ShiftDeviation(2);// lMS
		double fMS = S.getfMS(rMS, lMS);// fMS
		double GMS = S.getGMS(MS1, MS2, MS3, rMS, lMS);
		FileWriter Log;
		try {
			Log = new FileWriter("data/Log.txt", true);
			Log.write(sequenceFileName + "  " +new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date()) + ", " + MS1
					+ ", " + MS2 + ", " + MS3 + ", " + MS0 + ", " + rMS + ", " + lMS + ", " + fMS + ", " + GMS + "\n");
			Log.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
//
//		CodePermutation P=new CodePermutation();
//		P.generateCodes();
//		CodeFinder C=new CodeFinder();
//		C.RunCodeFinder(15);
//		new CodeEvaluation(P.calculateValues()).countBetterCodes();
		
	
	// write a statistic files

	// write statistic file 'Triplet Complete Transition'
	// try {
	// FileWriter fw = new FileWriter(new
	// File("data/TripletCompleteTransitionStat2012.csv"), true);
	// fw.write("\r\n Triplet Transiton Results calculated on "
	// + new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date()) + ":" +
	// "\r\n");
	// fw.write("Loaded Sequence filename: " + MainClass.sequenceFileName +
	// "\r\n" + "leading sequence" + ";"
	// + " sequence");
	//
	// String[] Bases = { "T", "C", "A", "G", "N" };
	// DecimalFormat df = new DecimalFormat("0.0000");
	// for (int x = 0; x < 4; x++) {
	// for (int y = 0; y < 4; y++) {
	// for (int z = 0; z < 4; z++) {
	//
	// for (int a = 0; a < 4; a++) {
	// for (int b = 0; b < 4; b++) {
	// for (int c = 0; c < 4; c++) {
	// fw.write("\r\n" + Bases[a] + Bases[b] + Bases[c] + ";" + Bases[x] +
	// Bases[y]
	// + Bases[z] + ";"
	// + df.format((double) tripletComleteTransition[x][y][z][a][b][c][0]) +
	// ";");
	// }
	// }
	// }
	//
	// }
	// }
	// }
	//
	// fw.write("\r\n \r\n");
	// fw.close();
	//
	// } catch (IOException e) {
	// System.out.println("Filewriter Error");
	// e.printStackTrace();
	// }

	// write statistic file 'Base apriori'
	// try {
	// FileWriter fw = new FileWriter(new File("data/BaseAprioriStat2012.csv"),
	// true);
	// fw.write("\r\n Base Apriori Results calculated on "
	// + new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date()) + ":" +
	// "\r\n");
	// fw.write("Loaded Sequence filename: " + sequenceFileName);
	// String[] Bases = { "T", "C", "A", "G", "N" };
	// DecimalFormat df = new DecimalFormat("0.0000");
	// for (int x = 0; x < 5; x++) {
	//
	// fw.write("\r\n" + Bases[x] + ";" + df.format(baseAprioriCount[x]) + ";");
	// }
	// fw.write("\r\n" + "Base Count: " + ";" + df.format(baseAprioriCount[5]));
	// fw.write("\r\n \r\n");
	// fw.close();
	//
	// } catch (IOException e) {
	// System.out.println("Filewriter Error");
	// e.printStackTrace();
	// }
	// // write statistic file 'Triplet apriori'
	// try {
	// FileWriter fw = new FileWriter(new
	// File("data/TripletAprioriStat2012.csv"), true);
	// fw.write("\r\n Triplet Apriori Results calculated on "
	// + new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date()) + ":" +
	// "\r\n");
	// fw.write("Loaded Sequence filename: " + sequenceFileName);
	//
	// String[] Bases = { "T", "C", "A", "G", "N" };
	// DecimalFormat df = new DecimalFormat("0.0000");
	// for (int x = 0; x < 5; x++) {
	// for (int y = 0; y < 5; y++) {
	// for (int z = 0; z < 5; z++) {
	// fw.write("\r\n" + Bases[x] + Bases[y] + Bases[z] + ";" +
	// df.format(tripletAprioriCount[x][y][z])
	// + ";");
	// }
	// }
	// }
	//
	// fw.write("\r\n \r\n");
	// fw.close();
	// } catch (IOException e) {
	// System.out.println("Filewriter Error");
	// e.printStackTrace();
	// }
	// // write statistic file 'Base Transition'
	// try {
	// FileWriter fw = new FileWriter(new
	// File("data/BaseTransitionStat2012.csv"), true);
	// fw.write("\r\n Base Transiton Results calculated on "
	// + new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date()) + ":" +
	// "\r\n");
	// fw.write("Loaded Sequence filename: " + sequenceFileName);
	//
	// String[] Bases = { "T", "C", "A", "G", "N" };
	// DecimalFormat df = new DecimalFormat("0.0000");
	// for (int x = 0; x < 5; x++) {
	// for (int y = 0; y < 5; y++) {
	// fw.write("\r\n" + Bases[x] + Bases[y] + ";" +
	// df.format(baseTransitionCount[x][y]) + ";");
	// }
	// }
	// fw.write("\r\n" + "Base Count: " + ";" + df.format(baseAprioriCount[5]));
	// fw.write("\r\n \r\n");
	// fw.close();
	//
	// } catch (IOException e) {
	// System.out.println("Filewriter Error");
	// e.printStackTrace();
	// }
	// // write statistic file 'Triplet Transition'
	// try {
	// FileWriter fw = new FileWriter(new
	// File("data/TripletTransitionStat2012.csv"), true);
	// fw.write("\r\n Triplet Transition Results calculated on "
	// + new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date()) + ":" +
	// "\r\n");
	// fw.write("Loaded Sequence filename: " + sequenceFileName + "\r\n" + ";" +
	// "leading T" + ";" + "following T"
	// + ";" + "leading C" + ";" + "following C" +";" + "leading A" + ";" +
	// "following A" +";" + "leading G" + ";"
	// + "following G" + ";" +"leading N" + ";" + "following N");
	//
	// String[] Bases = { "T", "C", "A", "G", "N" };
	// DecimalFormat df = new DecimalFormat("0.0000");
	// for (int x = 0; x < 5; x++) {
	// for (int y = 0; y < 5; y++) {
	// for (int z = 0; z < 5; z++) {
	// fw.write("\r\n" + Bases[x] + Bases[y] + Bases[z] + ";"
	// + df.format(tripletTransitionCount[x][y][z][0][0]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][0][1]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][1][0]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][1][1]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][2][0]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][2][1]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][3][0]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][3][1]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][4][0]) + ";"
	// + df.format(tripletTransitionCount[x][y][z][4][1]));
	// }
	// }
	// }
	// fw.write("\r\n \r\n");
	// fw.close();
	// } catch (IOException e) {
	// System.out.println("Filewriter Error");
	// e.printStackTrace();
	// }

	// 55417888 = 33MB Mixed
	// 51847843 = 158 MB Mixed Chromosom 7
	// 568815597 = 250 MB Mixed Chromosom 1
	// 671162122 = 32 MB Drosophila melanogaster chromosome 2R
	// Get One mixed String from Multi-Sequence File
	// List<DNASequence> Mixed=conn.LoadMixedFile();
	// StringBuilder builder = new StringBuilder();
	// for(DNASequence seq : Mixed) {
	// builder.append(seq.getSequenceAsString());
	// }
	//
	// String MixedSeq=builder.toString();
	// System.out.println("Mixed Length: "+MixedSeq.length());

	// // Berits try following
	// //Matrix and apriori!
	//
	// SequenceStats Stat=new SequenceStats(MixedSeq);
	//
	// baseAprioriWeights=Stat.getBase_aPriori();
	// tripletAprioriWeights=Stat.getTriplet_aPriori();
	// baseTransitionWeights=Stat.getBaseTransition();
	// tripletTransitionWeights=Stat.getTripletTransition();
	//
	// String[]Bases={"T","C","A","G"};
	// DecimalFormat df = new DecimalFormat("0.0000");
	// for (int x=0;x<4;x++){
	// for (int y=0;y<4;y++){
	// for (int z=0;z<4;z++){
	// System.out.println(Bases[x]+Bases[y]+Bases[z]+";"+df.format(tripletTransitionWeights[x][y][z][0][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][0][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][1][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][1][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][2][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][2][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][3][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][3][1]/4));//Warum
	// /4 ?
	// }
	// }
	// }
	//
	// // get MS1, MS0...
	//
	// TransitionTransversionBias=2;
	// setWeightings(true, true, true, true);
	// GeneCode G=new GeneCode(Code);
	// StabilityCalculator S=new StabilityCalculator(G);
	// double MS1=S.get_BaseDeviation(1);//MS1
	// double MS2=S.get_BaseDeviation(2);//MS2
	// double MS3=S.get_BaseDeviation(3);//MS3
	// double MS0=S.getMS0(MS1, MS2, MS3);//MS0
	// double rMS=S.get_ShiftDeviation(1);//rMS
	// double lMS=S.get_ShiftDeviation(2);//lMS
	// double fMS=S.getfMS(rMS, lMS);//fMS
	// double GMS=S.getGMS(MS1, MS2, MS3, rMS, lMS);
	// FileWriter Log;
	// try {
	// Log = new FileWriter("data/Log.txt", true);
	// Log.write("\r\n"+"Mixed Sequence: "+MS1+", "+MS2+", "+MS3+", "+MS0+",
	// "+rMS+", "+lMS+", "+fMS+", "+GMS+"\n");
	// Log.close();
	// } catch (IOException e) {
	// e.printStackTrace();
	// }
	//
	//
	// looking for better codes

	// TransitionTransversionBias=2;
	// CodePermutation P=new CodePermutation();
	// P.generateCodes();
	// setWeightings(false, true, false, true);
	// CodeFinder C=new CodeFinder();
	// C.RunCodeFinder(15);
	// new CodeEvaluation(P.calculateValues()).countBetterCodes();
	// Berits try end

	// double[] w=stat.getNucleotideDistribution(MixedSeq);
	// double[] tempfactors={w[0]/0.25,w[1]/0.25,w[2]/0.25,w[3]/0.25};
	// factors=tempfactors;
	// tweights=stat.getTripletDistribution(MixedSeq);

	// Batch Calculation to compare different sequences and their impact on
	// the Scores
	// List<String> Sequences=new ArrayList<String>();
	// Sequences.add("568815597"); //H.S 1
	// Sequences.add(568815596);//2
	// Sequences.add(568815595);//3
	// Sequences.add(568815594);//4
	// Sequences.add(568815593);//5
	// Sequences.add(568815592);//6
	// Sequences.add(568815591);//7
	// Sequences.add(568815590);//9
	// Sequences.add(568815588);//10
	// Sequences.add(568815587);//11
	// Sequences.add(568815586);//12
	// Sequences.add(568815585);//13
	// Sequences.add(568815584);//14
	// Sequences.add(568815583);//15
	// Sequences.add(568815582);//16
	// Sequences.add(568815581);//17
	// Sequences.add(568815580);//18
	// Sequences.add(568815579);//19
	// Sequences.add(568815578);//20
	// Sequences.add(568815577);//21
	// Sequences.add(568815576);//22
	// Sequences.add(568815575);//X
	// Sequences.add(568815574);//Y
	// Sequences.add(671162317);//D. Melanogaster 3L
	// Sequences.add(671162315);//D. Melanogaster 2R
	// Sequences.add(56384585);// E.Coli
	// Sequences.add(1537050);// HIV
	// Sequences.add("NC_001806"); //Herpes
	// Sequences.add("NC_001591"); //Papillomavirus
	// Sequences.add("NC_003461"); //Parainfluenza
	// Sequences.add("NC_001959"); //Norovirus G1
	// Sequences.add("NC_002031");//Gelbfieber
	// Sequences.add("NC_012532");//Zika
	// Sequences.add("NC_001477");//Dengue
	// Sequences.add("NC_004102");//Hepatitis C
	// Sequences.add("NC_002549");//Zaire Ebola

	// Sequences.add("mart_CDS_1_2_sequence"); // Berit writing...

	// for (String Str:Sequences){
	// DNASequence Seq1=conn.LoadFastaFile(Str);
	// SequenceStats Stat=new SequenceStats(Seq1.getSequenceAsString());
	// baseAprioriWeights=Stat.getBase_aPriori();
	// tripletAprioriWeights=Stat.getTriplet_aPriori();
	// baseTransitionWeights=Stat.getBaseTransition();
	// tripletTransitionWeights=Stat.getTripletTransition();
	// TransitionTransversionBias=2;
	// setWeightings(true, true, true, true);
	// GeneCode G=new GeneCode(Code);
	// StabilityCalculator S=new StabilityCalculator(G);
	// double MS1=S.get_BaseDeviation(1);//MS1
	// double MS2=S.get_BaseDeviation(2);//MS2
	// double MS3=S.get_BaseDeviation(3);//MS3
	// double MS0=S.getMS0(MS1, MS2, MS3);//MS0
	// double rMS=S.get_ShiftDeviation(1);//rMS
	// double lMS=S.get_ShiftDeviation(2);//lMS
	// double fMS=S.getfMS(rMS, lMS);//fMS
	// double GMS=S.getGMS(MS1, MS2, MS3, rMS, lMS);
	// FileWriter Log;
	// try {
	// Log = new FileWriter("data/Log.txt", true);
	// Log.write(Str+", "+MS1+", "+MS2+", "+MS3+", "+MS0+", "+rMS+",
	// "+lMS+", "+fMS+", "+GMS+"\n");
	// Log.close();
	// } catch (IOException e) {
	// e.printStackTrace();
	// }
	// }
	//
	//
	// DNASequence Seq1=conn.LoadFastaFile("568815597"); //chomosom 1
	// DNASequence Seq1=conn.LoadFastaFile("42"); // für ges. angeglichen
	// mit query key

	// DNASequence Seq1=conn.LoadFastaFile("mart_CDS_1_2_3sequence");
	// //Berit
	// SequenceStats Stat=new SequenceStats(Seq1.getSequenceAsString());
	//
	// baseAprioriWeights=Stat.getBase_aPriori();
	// tripletAprioriWeights=Stat.getTriplet_aPriori();
	// baseTransitionWeights=Stat.getBaseTransition();
	// tripletTransitionWeights=Stat.getTripletTransition();
	//
	// String[]Bases={"T","C","A","G"};
	// DecimalFormat df = new DecimalFormat("0.0000");
	// for (int x=0;x<4;x++){
	// for (int y=0;y<4;y++){
	// for (int z=0;z<4;z++){
	// System.out.println(Bases[x]+Bases[y]+Bases[z]+";"+df.format(tripletTransitionWeights[x][y][z][0][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][0][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][1][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][1][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][2][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][2][1]/4)+";"+df.format(tripletTransitionWeights[x][y][z][3][0]/4)+";"+df.format(tripletTransitionWeights[x][y][z][3][1]/4));//Warum
	// /4 ?
	// }
	// }
	// }
	// TransitionTransversionBias=2;
	//// CodePermutation P=new CodePermutation();
	//// P.generateCodes();
	// setWeightings(false, true, false, true);
	// CodeFinder C=new CodeFinder();
	// C.RunCodeFinder(15);
	// new CodeEvaluation(P.calculateValues()).countBetterCodes();
	// CodeFinder C=new CodeFinder();
	// setWeightings(true,true,true,true);
	// CodePermutation P=new CodePermutation();
	// new CodeEvaluation(P.calculateValues()).countBetterCodes();
	// C.RunCodeFinder(15);
	//
	// CodePermutation P=new CodePermutation();
	// P.generateCodes();
	//
	// CodeFinder C=new CodeFinder();
	// C.RunCodeFinder(20);
	// //Ohne Gewichtung
	// CodePermutation P=new CodePermutation();
	// new CodeEvaluation(P.calculateValues()).countBetterCodes();
	//

	// ###################################################################
	// ###################################################################
	// ###################################################################

	// Berit Runtime Test
	// Get the Java runtime
	Runtime runtime = Runtime.getRuntime();
	// Run the garbage collector
	runtime.gc();
	// Calculate the used memory
	long memory = runtime.totalMemory() - runtime
			.freeMemory();
	System.out.println("Used memory is bytes: "+memory);
	System.out.println("Used memory is megabytes: " +bytesToMegabytes(memory));
		// Runtime ende

	}

	// Set Weightings here to enable or disable them globally
	private static void setWeightings(boolean ba, boolean ta, boolean bt, boolean tt) {
		baseAprioriEnabled = ba;
		tripletAprioriEnabled = ta;
		baseTransitionEnabled = bt;
		tripletTransitionEnabled = tt;
	}

}
