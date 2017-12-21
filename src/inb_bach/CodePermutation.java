package inb_bach;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

import Objects.GeneCode;

public class CodePermutation {
	String[] Code={"Leu","Pro","His","Gln","Arg","Ile","Met","Thr","Asn","Lys","Ser","Val","Ala","Asp","Glu","Gly","Phe","Tyr","Cys","Trp"}; //natural code
	String[] RegC={"Phe","Leu","Ile","Met","Val","Ser","Pro","Thr","Ala","Tyr","His","Gln","Asn","Lys","Asp","Glu","Cys","Trp","Arg","Gly"};//Was ist das?
	int[] DiffCode={16,   -1,   3,    3,     7,    5,   -5,    0,    4,    8,    -8, -8,    -4,    -4,   -1,   -1,   2,    2,    -14,   -4}; //???
	private int CodeCount;
	private double[][] ValueBuffer;
	private int[] progress;
	private List<String> CodeBuffer;
	private int nextValue;
	private FileWriter codes;
	private FileWriter values;
	private Random rnd;
	private int Threads=10;
	private int ThreadsFinished;
	private class ThreadedCalculator extends Thread {
	    public void run() {
	    	GeneCode g=new GeneCode();
	    	StabilityCalculator S=new StabilityCalculator(g);
	    	while (true){
	    		int currentCode=getNextValue();
	    		
	    		//if (currentCode%10000==0)System.out.println("Generating Code "+currentCode);
	    		if(currentCode>=CodeBuffer.size())break;
	    		String[] rCode=CodeBuffer.get(currentCode).split(" ")[1].split("~");
	    		g.changeCode(rCode);
	    		//writes the dataset columns
	    		ValueBuffer[currentCode][0]=S.get_BaseDeviation(1);//MS1
	    		ValueBuffer[currentCode][1]=S.get_BaseDeviation(2);//MS2
	    		ValueBuffer[currentCode][2]=S.get_BaseDeviation(3);//MS3
	    		ValueBuffer[currentCode][3]=S.getMS0(ValueBuffer[currentCode][0], ValueBuffer[currentCode][1], ValueBuffer[currentCode][2]);//MS0
	    		ValueBuffer[currentCode][4]=S.get_ShiftDeviation(1);//rMS
	    		ValueBuffer[currentCode][5]=S.get_ShiftDeviation(2);//lMS
	    		ValueBuffer[currentCode][6]=S.getfMS(ValueBuffer[currentCode][4], ValueBuffer[currentCode][5]);//fMS
	    		ValueBuffer[currentCode][7]=S.getGMS(ValueBuffer[currentCode][0], ValueBuffer[currentCode][1], ValueBuffer[currentCode][2], ValueBuffer[currentCode][4], ValueBuffer[currentCode][5]);
	    		//Marks that this dataset if fully computed
	    		progress[currentCode]=1;
	    	}
	    	//System.out.println("Threads finished: "+(++ThreadsFinished));
	    }
	}
	
	private String[] ConvertCode(String[] Import){
		String[]C=new String[20];
		for (int i=0;i<20;i++){
			C[i+DiffCode[i]]=Import[i];
		}
		return C;
	}
	//Method to read exported Matlab codes from R.Geyer
	public void importCodes(){
		BufferedReader br;
		try {
			codes=new FileWriter("data/codes.txt");
			br = new BufferedReader(new FileReader("data/CodeImport.txt"));
			String line = null;
			CodeCount=0;
	    	while ((line = br.readLine()) != null) {
		    	String[]Temp=line.split(",");
		    	String[]ImportedCode=new String[20];
		    	for (int i=0;i<20;i++){
		    		ImportedCode[i]=RegC[Integer.parseInt(Temp[i])-1];
		    	}
		    	writeCodeLine(ConvertCode(ImportedCode));
	    	}
	    	br.close();
	    	codes.close();
		} catch (Exception e) {
			System.out.println("FEHLER:");
			e.printStackTrace();
			return;
		}
	}

	public double[][] calculateValues(){
		CodeCount=0;
		nextValue=0;
		ThreadsFinished=0;
		//Initialize Writer for Value Output
		CodeBuffer=new ArrayList<String>();
		try {
			values=new FileWriter("data/codeValues.csv");
		} catch (IOException e) {
			System.out.println("FileWriter Error!");
			return null;
		}
		System.out.println("Buffering Codes...");
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
			return null;
		}
		//Intitalize output Buffer Array
		ValueBuffer=new double[CodeBuffer.size()][8];
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
		return ValueBuffer;
	}
	private synchronized int getNextValue(){
		return nextValue++;
	}
	
	
	//Prints the values to file
	private void writeValueLine(int Number,double[] Line){
		//Writes Value Line To File
		try {
			values.write(Number+" "+Arrays.toString(Line)+'\n');
		} catch (IOException e) {
			System.out.println("FileWriter Error on "+Line);
		}
		if (Number%50000==0)System.out.println("Calculated "+Number+" Codes with stabilities");
	}
	
	//Writes random generated Code line to File
	private void writeCodeLine(String[] rCode){
		CodeCount++;{
		try {
			codes.write(CodeCount+" "+rCode[0]+"~"+rCode[1]+"~"+rCode[2]+"~"+rCode[3]+"~"+rCode[4]+"~"+rCode[5]+"~"+rCode[6]+"~"+rCode[7]+"~"+rCode[8]+"~"+rCode[9]+"~"+rCode[10]+"~"+rCode[11]+"~"+rCode[12]+"~"+rCode[13]+"~"+rCode[14]+"~"+rCode[15]+"~"+rCode[16]+"~"+rCode[17]+"~"+rCode[18]+"~"+rCode[19]+'\n');;
		} catch (IOException e) {
			System.out.println("FileWriter Error on "+CodeCount);
		}
		if (CodeCount%10000==0)System.out.println("Generated "+CodeCount+" Codes");
		return;
		}
	}
	
	
	//returns a random generated Code
	private String[] getRandomCode(){
		String[] rCode=new String[20];
		for (int i=0;i<20;i++){
			boolean found=false;
			int acid=0;
			while (found==false){
				acid=rnd.nextInt(20);
				if (!(Contains(rCode,Code[acid],i)))found=true;;
			}
			rCode[i]=Code[acid];
		}
		return rCode;
	}
	private boolean Contains(String[]Snip, String Code, int Level){
		for (int i=0;i<Level;i++){
			if (Snip[i].equals(Code))return true;
		}
		return false;
	}
	
	public void generateCodes(){
		rnd=new Random();
		File f=new File("data/codes.txt");
		if (f.exists()){
			System.out.println("codes.txt existiert bereits.");
			return;
		}
		try {
			codes=new FileWriter("data/codes.txt");
			GeneCode g=new GeneCode();
			//1st is Natural code
	    	writeCodeLine(Code);
	    	
	    	for (int i=0;i<1000000;i++){
		    	String[] rCode=getRandomCode();
		    	writeCodeLine(rCode);
		    }
	    	codes.close();
		} catch (IOException e) {
			System.out.println("FileWriter Error!");
			return;
		}
	}
}

