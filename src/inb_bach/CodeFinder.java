package inb_bach;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ThreadLocalRandom;

import Objects.GeneCode;

public class CodeFinder {
	String[] Code={"Leu","Pro","His","Gln","Arg","Ile","Met","Thr","Asn","Lys","Ser","Val","Ala","Asp","Glu","Gly","Phe","Tyr","Cys","Trp"};
	private int CodeCount;
	private double[]ValueBuffer;
	private int[] progress;
	private ArrayList<String> CodeBuffer;
	private int nextValue;
	private FileWriter codes;
	private FileWriter values;
	private int Threads=10;
	private int ThreadsFinished;
	private HashMap<String[], Double> GMSValues;
	private List<String[]> mutatedCodes;
	private List<String[]> bestCodes;
	
	//runs the code finder algorithm.
	public void RunCodeFinder(int MaxRecursion){
		for (int i=0;i<MaxRecursion;i++){
			calculateValues();
			printStatistics(i);
			if(i==MaxRecursion-1)break;
			mutateCodes(1000);
			clearDuplicates();
		}
	}
	//Removes duplicates from the set
	private void clearDuplicates(){
		LinkedHashSet<String> duplicateFree;
		List<String> RawStrings=new ArrayList<String>();
		for (String[] rCode:mutatedCodes){
			RawStrings.add(rCode[0]+"~"+rCode[1]+"~"+rCode[2]+"~"+rCode[3]+"~"+rCode[4]+"~"+rCode[5]+"~"+rCode[6]+"~"+rCode[7]+"~"+rCode[8]+"~"+rCode[9]+"~"+rCode[10]+"~"+rCode[11]+"~"+rCode[12]+"~"+rCode[13]+"~"+rCode[14]+"~"+rCode[15]+"~"+rCode[16]+"~"+rCode[17]+"~"+rCode[18]+"~"+rCode[19]);
		}
		duplicateFree=new LinkedHashSet<String>(RawStrings);
		System.out.println("Number of Duplicates cleared:"+(RawStrings.size()-duplicateFree.size()));
		System.out.println("Remaining Codes for next Run:"+duplicateFree.size());
		CodeCount=0;
		try {
			codes=new FileWriter("data/codes.txt");
			for (String S:duplicateFree){
				writeCodeLine(S);
			}
			codes.close();
			Thread.sleep(1000L);
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}
	}
	//generates permutations of the n best codes
	private void mutateCodes(int amount){
		bestCodes=new ArrayList<String[]>();
		//Extract best codes
		//System.out.println("First Key"+GMSValues.firstKey());
		System.out.println("Extracting best Codes...");
		List<Double> ValuesTemp=new ArrayList<Double>(GMSValues.values());
		Collections.sort(ValuesTemp);
		double threshold;
		if (amount>=ValuesTemp.size()){
			threshold=ValuesTemp.get(ValuesTemp.size()-1);
		}else{
			threshold=ValuesTemp.get(amount);
		}
		for (Entry<String[], Double> e:GMSValues.entrySet()){
			if (e.getValue()<=threshold){
				bestCodes.add(e.getKey());
			}
		}
		//Permutation
		System.out.println("Permutating the best codes to find even better ones...");
		mutatedCodes=new ArrayList<String[]>();
		mutatedCodes.addAll(bestCodes);
		for (int i=0;i<bestCodes.size();i++){
			String[] Original=bestCodes.get(i);
			for (int pos1=0;pos1<20;pos1++){
				for (int pos2=pos1+1;pos2<20;pos2++){
					String[] NewCode=new String[20];
					System.arraycopy(Original, 0, NewCode, 0, 20);
					NewCode[pos1]=Original[pos2];
					NewCode[pos2]=Original[pos1];
					mutatedCodes.add(NewCode);
				}
			}
		}
		//Also include 5-Switch Permutations of the best codes to avoid local minimum
		for (int i=0;i<bestCodes.size();i++){
			String[] Original=bestCodes.get(i);
			String[] Changed=new String[20];
			System.arraycopy(Original, 0, Changed, 0, 20);
			int x=0;
			while (x<5){
				Random ran = new Random();
				int x1 = ran.nextInt(20);
				int x2 = ran.nextInt(20);
				if (x1!=x2){
					x++;
					String[] NewCode=new String[20];
					System.arraycopy(Changed, 0, NewCode, 0, 20);
					NewCode[x1]=Changed[x2];
					NewCode[x2]=Changed[x1];
					Changed=NewCode;
				}
			}
			mutatedCodes.add(Changed);
		}
	}
	//Prints the statistics of the code set on the console
	private void printStatistics(int i){
		double mean=0;
		double min=999999999;
		double max=-999999999;
		int cnt=0;
		for (String[] rCode:GMSValues.keySet()){
			double x=GMSValues.get(rCode);
			if (x==0)System.out.println("NULL");
			mean=mean+x;
			if (x<min)min=x;
			if (x>max)max=x;
			cnt++;
		}
		mean=mean/(double)cnt;
		System.out.println("Minimum: "+min);
		System.out.println("Maximum: "+max);
		System.out.println("Mean Value: "+mean);
		try {
	    	 FileWriter fw = new FileWriter(new File("data/GreedyResults.log"), true);
	    	 if (i==0)fw.write("Date, Run,  Minimum,  Maximum,  Mean "+"\n");
	    	 fw.write(new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date())+", "+i+", ");
	    	 fw.write(min+", "+max+", "+mean+"\n");
	    	 fw.close();
	    } catch (IOException e) {
			System.out.println("Filewriter Error");
			e.printStackTrace();
		}
		if(mean<min){
			System.out.println("FEHLER: Mean too small");
			System.exit(0);
		}
		
	}
	
	private synchronized int getNextValue(){
		return nextValue++;
	}
	private class ThreadedCalculator extends Thread {
	    public void run() {
	    	GeneCode g=new GeneCode();
	    	StabilityCalculator S=new StabilityCalculator(g);
	    	while (true){
	    		int currentCode=getNextValue();
	    		if(currentCode>=CodeBuffer.size())break;
	    		String[] rCode=CodeBuffer.get(currentCode).split(" ")[1].split("~");
	    		g.changeCode(rCode);
	    		double GMS=S.getGMS(S.get_BaseDeviation(1), S.get_BaseDeviation(2), S.get_BaseDeviation(3), S.get_ShiftDeviation(1), S.get_ShiftDeviation(2));
	    		GMSValues.put(rCode, GMS);
	    		ValueBuffer[currentCode]=GMS;
	    		//Marks that this dataset if fully computed
	    		progress[currentCode]=1;
	    	}
	    }
	}
	//Calculates GMS scores for each code in the set.
	private boolean calculateValues(){
		CodeCount=0;
		nextValue=0;
		ThreadsFinished=0;
		CodeBuffer=null;
		ValueBuffer=null;
		progress=null;
		GMSValues=new HashMap<String[], Double>();
		//Initialize Writer for Value Output
		try {
			values=new FileWriter("data/codeValues.csv");
			//Cache all Codes into Buffer
			CodeBuffer=new ArrayList<String>();
			
		} catch (IOException e) {
			System.out.println("FileWriter Error!");
			return false;
		}
		System.out.println("Buffering Codes..."); //TODO hier bricht das Programm ab! 
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader("data/codes.txt"));
			String line = null;
	    	while ((line = br.readLine()) != null) {
	    		CodeBuffer.add(line);
	    	}
	    	br.close();
		} catch (Exception e) {
			System.out.println("FEHLER: Filereader Error");
			e.printStackTrace();
			return false;
		}
		//Intitalize output Buffer Array
		ValueBuffer=new double[CodeBuffer.size()];
		progress=new int[CodeBuffer.size()];
		System.out.println("Start calculation with " +Threads+" Threads...");
		//Start the Calculation Threads
		for (int i=0;i<Threads;i++){
			new ThreadedCalculator().start();
			try {
				Thread.sleep(1);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		//Write Task to save the Results
		//Sleeps for 10ms if it is too fast
		for (int i=0;i<CodeBuffer.size();i++){
			int timeout=0;
			boolean repeat=false;
			if (progress[i]==0){
				if (ThreadsFinished==Threads)break;
				try {
					Thread.sleep(100);
					if (progress[i]!=0)repeat=true;;
					timeout+=100;
				} catch (InterruptedException e) {}
				if (timeout>3000)System.out.println("Timeout: "+timeout);
			}
			if (repeat){
				i--;
			}else{
				writeValueLine(i+1, ValueBuffer[i]);
			}
		}
		try {
			values.close();
		} catch (IOException e) {}
		return true;
	}
	private void writeValueLine(int Number,double Line){
		//Writes Value Line To File
		try {
			values.write(Number+" "+Line+'\n');
		} catch (IOException e) {
			System.out.println("FileWriter Error on "+Line);
		}
		if (Number%50000==0)System.out.println("Calculated "+Number+" Codes with stabilities");
	}
	
	//Writes random generated Code line to File
	private void writeCodeLine(String rCode){
		CodeCount++;{
		try {
			codes.write(CodeCount+" "+rCode+'\n');;
		} catch (IOException e) {
			System.out.println("FileWriter Error on "+CodeCount);
		}
		return;
		}
	}
}
