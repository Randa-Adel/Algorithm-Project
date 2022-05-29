#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<list>
#include<algorithm>
using namespace std;
vector<string> String_Suffix;
vector<string> Reads;
vector<string> ID;
vector<long long> Sorted_Index_Of_Reads;
vector<string> choromes;
struct SuffixArray
{
	string Seq;
	int index;
};
vector<SuffixArray>SArray;
void MakeSuffixArray(string seq)
{
	SuffixArray s;
	SArray.push_back(s);
	int k = 0;
	seq += "$";
	long long Last_index = seq.find("$");
	while (k != Last_index)
	{
		string temp = "";
		for (long long i = k; i < Last_index; i++)
		{
			temp += seq[i];
		}
		temp += "$";
		SArray.push_back(s);
		SArray[k].Seq = temp;
		SArray[k].index = k;
		k++;
	}

	SArray[k].Seq = "$";
	SArray[k].index = k;
	/*بنقل هنا السيكونس من لاسفكس اراي الكبير واحطها في فيكتور من الاسترينج للسفكس ارراي*/
	for (long long i = 0; i < SArray.size(); i++)
	{
		String_Suffix.push_back(SArray[i].Seq);
	}
}
void SortingSuffixArray()
{
	for (int i = 0; i < SArray.size(); i++)
	{
		for (int j = 0; j < SArray.size() - 1; j++)
		{
			if (SArray[i].Seq < SArray[j].Seq)
			{
				string temp = SArray[i].Seq;
				SArray[i].Seq = SArray[j].Seq;
				SArray[j].Seq = temp;
				int temp2 = SArray[i].index;
				SArray[i].index = SArray[j].index;
				SArray[j].index = temp2;
			}
		}
	}

}
void StoreSuffix()
{
	string line;
	ifstream file;
	file.open("SuffixArray.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!"<<endl;
		while (!file.eof())
		{
			getline(file, line);
			String_Suffix.push_back(line);
			
		}
	}
	file.close();
}
void StoreReads()
{
	string line;
	ifstream file;
	file.open("Reads.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			Reads.push_back(line);

		}
	}
	file.close();
}
void StoreID()
{
	string line;
	ifstream file;
	file.open("ID.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			ID.push_back(line);

		}
	}
	file.close();
}
void searchAboutRead()
{
	vector<string> TReads;
	for (size_t i = 0; i < Reads.size(); i++)
	{
		TReads.push_back(Reads[i] + "$");
	}
	for (int i = 0; i < Reads.size() ;i++)
	{
		bool flag = false;
		string read = TReads[i];
		int low = 0;
		int high = String_Suffix.size();
		while (low <= high) {
			int mid = low + (high - low) / 2;
			if (String_Suffix[mid].compare(read) == 0)
			{
				Sorted_Index_Of_Reads.push_back(mid);
				flag = true;				
			}
			if (String_Suffix[mid].compare(read)<0)
				low = mid + 1;
			else
				high = mid - 1;
		}	
		read = Reads[i];
		if (flag == false)
		{
			for (long long i = 0; i < String_Suffix.size(); i++)
			{
				if (String_Suffix[i].find(read) != -1)
				{
					Sorted_Index_Of_Reads.push_back(String_Suffix[i].find(read));
					flag = true;
					break;
				}
			}			
		}
		if (flag == false)
		{
			Sorted_Index_Of_Reads.push_back(-1);
		}
	}
	
}
void MakeSAMFormat(string seq)
{
	ofstream out("SAMFile.txt");
	out << "READ \t\t\t\t\t\t\t\t\t\t\t" << "READ-ID \t\t\t\t\t\t\t\t\t\t\t" << "Chromosome Number\t\t\t\t" << "Position In Refrence Genome " << endl;
	for (long long i = 0; i < Reads.size(); i++)
	{
		/*if (choromes[0].find(Reads[i])!=-1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 0\t\t\t\t\t" << seq.find(Reads[i]) << endl;
		}
		else if (choromes[1].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 1\t\t\t\t\t" << seq.find(Reads[i]) << endl;
		}
		else if (choromes[2].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 2\t\t\t\t\t" << seq.find(Reads[i]) << endl;
		}
		else if (choromes[3].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 3\t\t\t\t\t" << seq.find(Reads[i]) << endl;
		}
		else if (choromes[4].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 4\t\t\t\t\t" << seq.find(Reads[i]) << endl;
		}
		else if (choromes[5].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 5\t\t\t\t\t" << seq.find(Reads[i]) << endl;
		}*/
		if (choromes[0].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 0\t\t\t\t\t" << seq.find(Reads[i])+Reads[i].length() << endl;
		}
		else if (choromes[1].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 1\t\t\t\t\t" << seq.find(Reads[i]) + Reads[i].length() << endl;
		}
		else if (choromes[2].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 2\t\t\t\t\t" << seq.find(Reads[i]) + Reads[i].length() << endl;
		}
		else if (choromes[3].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 3\t\t\t\t\t" << seq.find(Reads[i]) + Reads[i].length() << endl;
		}
		else if (choromes[4].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 4\t\t\t\t\t" << seq.find(Reads[i]) + Reads[i].length() << endl;
		}
		else if (choromes[5].find(Reads[i]) != -1)
		{
			out << Reads[i] << "\t\t\t" << ID[i] << "\t\t\t\t\t" << "\t\t 5\t\t\t\t\t" << seq.find(Reads[i]) + Reads[i].length() << endl;
		}
		
	}
	out.close();


}
string CollectChromosomes()
{

	string line;
	ifstream file;
	string seq = "";
	file.open("chr01.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!"<<endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
			seq += line;
			}
		}
	}
	choromes.push_back(seq);
	file.close();
	file.open("chr02.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	choromes.push_back(seq);
	file.close();
	file.open("chr03.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	choromes.push_back(seq);
	file.close();
	file.open("chr04.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	choromes.push_back(seq);
	file.close();
	file.open("chr05.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	choromes.push_back(seq);
	file.close();
	file.open("chr06.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	choromes.push_back(seq);
	file.close();
	/*file.open("chr06.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	file.close();
	file.open("chr07.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	file.close();
	file.open("chr08.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	file.close();
	file.open("chr09.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}
	file.close();
	file.open("chr10.txt", ios::out | ios::in);
	if (!file) {
		cout << "File not acessible!";
	}
	else {
		cout << "File reached successfully!" << endl;
		while (!file.eof())
		{
			getline(file, line);
			if (line[0] == '>') { continue; }
			else
			{
				seq += line;
			}
		}
	}*/
	//file.close();
	
	ofstream out("Refrence Genome.txt");
	out << seq;
	out.close();

	return seq;
}
int main()
{
    string seq = CollectChromosomes();
	MakeSuffixArray(seq);
	SortingSuffixArray();
	StoreReads();
	StoreID();
	searchAboutRead();
	ofstream out("ReadsInSuffixArray.txt");
	out << "READ \t\t\t\t\t\t\t\t\t\t\t" << "READ-ID \t\t\t\t\t\t\t\t\t\t\t" << "\t\t\t\tINDEX IN SUFFIX ARRAY \t\t\t\t" << endl;
	for (long long i = 0; i < Sorted_Index_Of_Reads.size(); i++)
	{
		out <<Reads[i]<<"\t\t\t\t\t"<<ID[i]<<"\t\t\t\t\t"<< Sorted_Index_Of_Reads[i] << endl;
		
	}
	out.close();
	MakeSAMFormat(seq);
	return 0;
}