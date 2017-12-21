package inb_bach;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

public class CodeEvaluation {
	private double[][] values;
	
	public CodeEvaluation(double[][] values){
		this.values=values;
	}
	//Counts how many of the code in the current set are better than the first code in the set.
	//To work properly, the natural code should be the first code.
	//Needs the generated values matrix.
	public int[] countBetterCodes(){
		System.out.println("Evaluating Results and counting better codes found...");
		int[] betterCodes=new int[values[0].length+1];
		for (int i=0;i<values.length;i++){
			boolean isCompleteBetter=true;
			//Natural code needs to be #1
			for (int x=0;x<values[0].length;x++){
				if(values[i][x]<values[0][x]){
					betterCodes[x]++;
				}else{
					isCompleteBetter=false;
				}
			}
			if(isCompleteBetter){
				betterCodes[betterCodes.length-1]++;
				System.out.println("Allover Better Code found: "+i);
			}
		}
		System.out.println("Number of better codes found:");
		System.out.println(Arrays.toString(betterCodes)); 
	    try {
	    	 FileWriter fw = new FileWriter(new File("data/EvaluationResults.log"), true);
	    	 fw.write("Results calculated on "+new SimpleDateFormat("dd.MM.yyyy HH:mm:ss").format(new Date())+":"+"\n");
	    	 fw.write("   ["+MainClass.baseAprioriEnabled+ ","+ MainClass.tripletAprioriEnabled+ ","+ MainClass.baseTransitionEnabled+ ","+ MainClass.tripletTransitionEnabled+","+MainClass.TransitionTransversionBias+"]"+"\n");
	    	 fw.write("   "+Arrays.toString(betterCodes)+"\n");
	    	 fw.write("\n");
	    	 fw.close();
	    } catch (IOException e) {
			System.out.println("Filewriter Error");
			e.printStackTrace();
		}
		return betterCodes;
	}
//Prüft ob der Wert der orginal Sequenz größer ist. Wenn ja: neuer code besser. speichert und zählt bessere Codes. Gibt sie auf der Kosole aus und speichert sie in data/EvaluationResults.log  
}
