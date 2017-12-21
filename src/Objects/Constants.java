package Objects;

public class Constants {
	//Holds data for Polar requirements
	public static double getPolarReq(String AminoAcid){
		switch(AminoAcid){
		case "Ala":
			return 7.0;
		case "Asp":
			return 13.0;
		case "Cys":
			return 4.8;
		case "Gln":
			return 8.6;
		case "His":
			return 8.4;
		case "Leu":
			return 4.9;
		case "Met":
			return 5.3;
		case "Pro":
			return 6.6;
		case "Thr":
			return 6.6;
		case "Tyr":
			return 5.4;
		case "Arg":
			return 9.1;
		case "Asn":
			return 10.0;
		case "Glu":
			return 12.5;
		case "Gly":
			return 7.9;
		case "Ile":
			return 4.9;
		case "Lys":
			return 10.1;
		case "Phe":
			return 5.0;
		case "Ser":
			return 7.5;
		case "Trp":
			return 5.2;
		case "Val":
			return 5.6;
		}
		System.out.println("Unrecognized Amino Acid: "+ AminoAcid);
		return 0;
	}
}
