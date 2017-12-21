package inb_bach;

import Objects.Constants;
import Objects.GeneCode;

public class StabilityCalculator {
	private GeneCode Code;
	private String[]Bases={"T","C","A","G"};
	
	private boolean baseAprioriWeighting=false;
	private double[] baseAprioriWeights;
	
	private boolean tripletaPrioriWeighting=false;
	private double[][][] tripletAprioriWeights;
	
	private boolean baseTransitionWeighting=false;
	private double[][] baseTransitionWeights;
	
	private boolean tripletTransitionWeighting=false;
	private double[][][][][]tripletTransitionWeights;
	
	private int Bias=1;
	/* Deviationmode:
	 * 1=MS1 - Default
	 * 2=MS2
	 * 3=MS3
	 */
	public StabilityCalculator(GeneCode Code){
		this.Code=Code;
		this.Bias=MainClass.TransitionTransversionBias;
		//Changes the a Priori Weights for single Bases
		if (MainClass.baseAprioriEnabled){
			baseAprioriWeighting=true;
			baseAprioriWeights=MainClass.baseAprioriWeights;
		}
		//Changes the a Priori Weights for Triplets
		if (MainClass.tripletAprioriEnabled){
			tripletaPrioriWeighting=true;
			tripletAprioriWeights=MainClass.tripletAprioriWeights;
		}
		//Changes the Base-Transition Matrix
		if (MainClass.baseTransitionEnabled){
			baseTransitionWeighting=true;
			baseTransitionWeights=MainClass.baseTransitionWeights;
		}
		//Changes the Triplet-Transition Matrix used for Shift Calculations
		if (MainClass.tripletTransitionEnabled){
			tripletTransitionWeighting=true;
			tripletTransitionWeights=MainClass.tripletTransitionWeights;
		}
	}
	//Changes the Code used for Calculations
	public void ChangeCode(String[] Mapping){
		Code.changeCode(Mapping);
	}
	//Returns Deviations calculated with all current set Parameters //Returns MS1, MS2 oder MS3
	//Only for single base Mutations
	//1=MS1
	//2=MS2
	//3=MS3
	public double get_BaseDeviation(int Modus){
		if(!(Modus>=1&&Modus<=3)){
			System.out.println("Bad Mode, Allowed: 1-3");
			return 0;
		}
		double deviation=0.0;
		for (int i=0;i<4;i++){
			String a=Bases[i];
			for (int j=0;j<4;j++){
				String b=Bases[j];
				for (int k=0;k<4;k++){
					String c=Bases[k];
					String Amino=Code.getAminoAcid(a+b+c);
					if (Amino.length()!=3)continue; //Filters Stop Codons
					double Polar1=Constants.getPolarReq(Amino);
					Double Diff=0.0;
					for (int m=0;m<4;m++){
						String x=Bases[m];
						String Amino2="";
						switch (Modus){
						case 1:
							Amino2=Code.getAminoAcid(x+b+c);
							break;
						case 2:
							Amino2=Code.getAminoAcid(a+x+c);
							break;
						case 3:
							Amino2=Code.getAminoAcid(a+b+x);
							break;
						}
						if (Amino2.length()!=3)continue; //Filters Stop Codons
						double Polar2=Constants.getPolarReq(Amino2);
						double difference=(Polar1-Polar2)*(Polar1-Polar2);
						if (baseAprioriWeighting){
							difference=difference*baseAprioriWeights[m];
						}
						if (tripletaPrioriWeighting){
							difference=difference*tripletAprioriWeights[i][j][k];
						}
						if(Bias!=1){
							String from="";
							switch (Modus){
							case 1:
								from=a;
								break;
							case 2:
								from=b;
								break;
							case 3:
								from=c;
								break;
							}
							if (isTransition(from, x)){
								difference=difference*Bias;
							}
						}
						Diff+=difference;
					}
					deviation=deviation+Diff;
				}
			}
		}
		switch (Modus){
		case 1:
			deviation=deviation/(174+((Bias-1)*58));
			break;
		case 2:
			deviation=deviation/(176+((Bias-1)*60));
			break;
		case 3:
			deviation=deviation/(176+((Bias-1)*60));
			break;
		}
		return deviation;
	}
	//returns deviation for Shift Mutations
	//1=Right (+1), 2=Left (-1)
	public double get_ShiftDeviation(int Modus){
		if(!(Modus>=1&&Modus<=2)){
			System.out.println("Bad Mode, Allowed: 1-2");
			return 0;
		}
		double deviation=0.0;
		for (int i=0;i<4;i++){
			String a=Bases[i];
			for (int j=0;j<4;j++){
				String b=Bases[j];
				for (int k=0;k<4;k++){
					String c=Bases[k];
					String Amino=Code.getAminoAcid(a+b+c);
					if (Amino.length()!=3)continue; //Filtert Stop Codons
					double Polar1=Constants.getPolarReq(Amino);
					Double Diff=0.0;
					for (int m=0;m<4;m++){
						String x=Bases[m];
						String Amino2="";
						switch (Modus){
						case 1:
							Amino2=Code.getAminoAcid(b+c+x);
							break;
						case 2:
							Amino2=Code.getAminoAcid(x+a+b);
							break;
						}
						if (Amino2.length()!=3)continue; //Filtert Stop Codons
						double Polar2=Constants.getPolarReq(Amino2);
						double difference=(Polar1-Polar2)*(Polar1-Polar2);
						if (baseAprioriWeighting){
							difference=difference*baseAprioriWeights[m];
						}
						if (tripletaPrioriWeighting){
							difference=difference*tripletAprioriWeights[i][j][k];
						}
						if (baseTransitionWeighting){
							switch(Modus){
							case 1:
								difference=difference*baseTransitionWeights[k][m];
								break;
							case 2:
								difference=difference*baseTransitionWeights[m][i];
								break;
							}
						}
						if(tripletTransitionWeighting){
							switch(Modus){
							case 1:
								difference=difference*tripletTransitionWeights[i][j][k][m][1];
								break;
							case 2:
								difference=difference*tripletTransitionWeights[i][j][k][m][0];
								break;
							}
						}
						Diff+=difference;
					}
					deviation=deviation+Diff;
				}
			}
		}
		//232 Possible Mutations...
		deviation=deviation/232;
		return deviation;
	}
	
	//returns true if the Mutation was a Transition
	private boolean isTransition(String from, String to){
		if(from.equalsIgnoreCase("A")&&to.equalsIgnoreCase("G"))return true;
		if(from.equalsIgnoreCase("G")&&to.equalsIgnoreCase("A"))return true;
		if(from.equalsIgnoreCase("C")&&to.equalsIgnoreCase("T"))return true;
		if(from.equalsIgnoreCase("T")&&to.equalsIgnoreCase("C"))return true;
		return false;
	}
	
	public double getWMS0(int Bias, double WMS1, double WMS2, double WMS3){
		WMS1=WMS1*(174+(Bias*58)); // warum *(174+(Bias*58) ???  es gilt: *(174+((Bias-1)*58) =F
		WMS2=WMS2*(176+(Bias*60));
		WMS3=WMS3*(176+(Bias*60));
		return (WMS1+WMS2+WMS3)/((174+(Bias*58))+176+(Bias*60)+176+(Bias*60)); // Das müsste doch: (WMS1+WMS2+WMS3)/((174+((Bias-1)*58))+176+((Bias-1)*60)+176+((Bias-1)*60))
	}
	
	public double getMS0(double MS1, double MS2, double MS3){
		MS1=MS1*174; // warum *F ???
		MS2=MS2*176;
		MS3=MS3*176;
		return(MS1+MS2+MS3)/(174+176+176);
	}
	
	public double getfMS(double rMS, double lMS){
		rMS=rMS*232; // warum *F ??? 
		lMS=lMS*232;
		return(rMS+lMS)/(232+232);
	}
	public double getGMS(double MS1, double MS2, double MS3, double rMS, double lMS){
		MS1=MS1*174; // Gleiche Frage wie oben
		MS2=MS2*176;
		MS3=MS3*176;
		rMS=rMS*232;
		lMS=lMS*232;
		return(MS1+MS2+MS3+rMS+lMS)/(174+176+176+232+232); // Bias(=T) wird hier nicht beachtet! 
	}
}
