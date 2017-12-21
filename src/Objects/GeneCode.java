package Objects;

public class GeneCode {
	//Each array holds all triplets coding for a specific amino acid.
	String[]Codes1={"CTT","CTC","CTA","CTG","TTA","TTG"};//Leu
	String[]Codes2={"CCT","CCC","CCA","CCG"};//Pro
	String[]Codes3={"CAT","CAC"};//His
	String[]Codes4={"CAA","CAG"};//Gln
	String[]Codes5={"CGT","CGC","CGA","CGG", "AGA", "AGG"};//Arg
	String[]Codes6={"ATT","ATC","ATA"};//Ile
	String[]Codes7={"ATG"};//Met & START
	String[]Codes8={"ACT","ACC","ACA","ACG"};//Thr
	String[]Codes9={"AAT","AAC"};//Asn
	String[]Codes10={"AAA","AAG"};//Lys
	String[]Codes11={"AGT","AGC","TCT","TCC","TCA","TCG"};//Ser
	String[]Codes12={"GTT","GTC","GTA","GTG"};//Val
	String[]Codes13={"GCT","GCC","GCA","GCG"};//Ala
	String[]Codes14={"GAT","GAC"};//Asp
	String[]Codes15={"GAA","GAG"};//Glu
	String[]Codes16={"GGT","GGC","GGA","GGG"};//Gly
	String[]Codes17={"TTT","TTC"};//Phe
	String[]Codes18={"TAT","TAC"};//Tyr
	String[]Codes19={"TGT","TGC"};//Cys
	String[]Codes20={"TGG"};//Trp
	//String[]Codes21={"TAA","TAG","TGA"};//Stops
	String[][]Allcodes={Codes1,Codes2,Codes3,Codes4,Codes5,Codes6,Codes7,Codes8,Codes9,Codes10,Codes11,Codes12,Codes13,Codes14,Codes15,Codes16,Codes17,Codes18,Codes19,Codes20};
	String[]Mapping;
	public GeneCode(){
		//Natural Code
		String[] Code={"Leu","Pro","His","Gln","Arg","Ile","Met","Thr","Asn","Lys","Ser","Val","Ala","Asp","Glu","Gly","Phe","Tyr","Cys","Trp"};
		Mapping=Code;
	}
	public GeneCode(String[]Mapping){
		//Can be used to instantiate any generated code
		this.Mapping=Mapping;
	}
	public void changeCode(String[] Mapping){
		//Can be used to instantiate any generated code
		this.Mapping=Mapping;
	}
	public String getAminoAcid(String Codon){
		if(Codon.length()!=3){
			System.out.println("Unrecognized Codon");
			return "ERROR";
		}
		for (int i=0;i<20;i++){
			for (int j=0;j<Allcodes[i].length;j++){
				if (Allcodes[i][j].equals(Codon)){
					return Mapping[i];
				}
			}
		}
		return "ERROR 2";
	}
}
