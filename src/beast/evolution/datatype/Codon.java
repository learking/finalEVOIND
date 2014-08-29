package beast.evolution.datatype;

import beast.evolution.datatype.DataType.Base;

public class Codon extends Base {

	public Codon () {
		stateCount = 61;
		codeLength = 3;
		codeMap = "AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTT";
		mapCodeToStateSet = new int[61][1];
        for (int i = 0; i < 61; i++) {
            mapCodeToStateSet[i][0] = i;
        }
	}

	@Override
	public String getDescription() {
		return "codon";
	}

}