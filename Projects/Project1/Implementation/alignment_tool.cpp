#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <limits>

using namespace std;

/*  parse_sequence_file :It takes a string as input, which is the name of  the file. 
The function returns a vector of characters, which are the contents of the file,
 but with all non-alphabetical characters removed.
*/

vector<char> parse_sequence_file(const string& file_path) {
    vector<char> sequences;
    string header, line;
    ifstream seqFile(file_path);
    if (!seqFile) {
        cerr << "Error in opening file, please enter a valid file name." << endl;
        return sequences;
    }
    getline(seqFile, header);
    getline(seqFile, line);
    for (char c : line) {
        if (isalpha(c)) {
            sequences.push_back(c);
        }
    }
    seqFile.close();
    return sequences;
}

map<pair<char, char>, int> parse_matrix_file(const string& file_path) {
    map<pair<char, char>, int> matrix;
    string line;
    char amino;
    size_t pos, row_counter;
    vector<char> amino_names;

    ifstream matrix_file(file_path);
    if (!matrix_file.is_open()) {
        cerr << "Error in opening file, please enter a valid file name." << endl;
        return matrix;
    }
    getline(matrix_file, line);
    getline(matrix_file, line);
    istringstream ss(line);

    while (ss >> amino) {
        if (ss.peek() == ',') {
            ss.ignore();
        }
        amino_names.push_back(amino);
    }
    row_counter = 0;
    while (getline(matrix_file, line)) {
        vector<string> values;
        pos = 0;
        while ((pos = line.find(',')) != string::npos) {
            values.push_back(line.substr(0, pos));
            line.erase(0, pos + 1);
        }
        values.push_back(line);
        for (size_t i = 0; i < values.size(); ++i) {
            matrix[{amino_names[row_counter], amino_names[i]}] = stoi(values[i]);
        }
        row_counter++;
    }
    return matrix;
}

// Implementation of the global_alignment

void global_alignment(const vector<char>& test1, const vector<char>& test2,
                      const map<pair<char, char>, int>& matrix, int gap_penalty) {
    // Initialize dimensions
    int m = test1.size();
    int n = test2.size();

    // Initialize a matrix to store the alignment scores
    vector<vector<int>> OPTMatrix(m + 1, vector<int>(n + 1, 0));

    // Initialize the first column and first row with gap penalties
    for (int i = 1; i <= m; ++i) {
        OPTMatrix[i][0] = OPTMatrix[i - 1][0] + gap_penalty;
    }

    for (int j = 1; j <= n; ++j) {
        OPTMatrix[0][j] = OPTMatrix[0][j - 1] + gap_penalty;
    }

    // Fill the matrix with scores
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match_score = OPTMatrix[i - 1][j - 1] + matrix.at(make_pair(test1[i - 1], test2[j - 1]));
            int gap1_score = OPTMatrix[i - 1][j] + gap_penalty;
            int gap2_score = OPTMatrix[i][j - 1] + gap_penalty;
            OPTMatrix[i][j] = max({match_score, gap1_score, gap2_score});
        }
    }

    // Traceback to find the aligned sequences
    vector<char> alignment_first_sequence, alignment_second_sequence;
    int i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && OPTMatrix[i][j] == OPTMatrix[i - 1][j] + gap_penalty) {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), test1[i - 1]);
            alignment_second_sequence.insert(alignment_second_sequence.begin(), '-');
            i -= 1;
        } else if (j > 0 && OPTMatrix[i][j] == OPTMatrix[i][j - 1] + gap_penalty) {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), '-');
            alignment_second_sequence.insert(alignment_second_sequence.begin(), test2[j - 1]);
            j -= 1;
        } else {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), test1[i - 1]);
            alignment_second_sequence.insert(alignment_second_sequence.begin(), test2[j - 1]);
            i -= 1;
            j -= 1;
        }
    }

    int total_score = OPTMatrix[m][n];

    // Print the aligned sequences
    cout << "Alignment of Sequence 1: ";
    for (const char& s : alignment_first_sequence) {
        cout << s;
    }
    cout << endl;

    cout << "Alignment of Sequence 2: ";
    for (const char& s : alignment_second_sequence) {
        cout << s;
    }
    cout << endl;

    // Print the OPT Matrix
    cout << "OPT Matrix: " << endl;
    for (const std::vector<int>& row : OPTMatrix) {
        for (int value : row) {
            std::cout << value << ", ";
        }
        std::cout << std::endl;
    }

    // Print optimal score
    cout << "Optimal alignment Score: " << total_score << endl;
}

// Implementation of the local_alignment

void local_alignment(const vector<char>& seq1, const vector<char>& seq2,
                     const map<pair<char, char>, int>& matrix, int gap_penalty) {
    int m = seq1.size();
    int n = seq2.size();

    // Initialize the table
    vector<vector<int>> OPTMatrix(m + 1, vector<int>(n + 1, 0));

    // Initialize variables to keep track of the maximum score and its position
    int max_score = 0;
    int max_i = 0, max_j = 0;

    // Fill in the table
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match_score = OPTMatrix[i - 1][j - 1] + matrix.at(make_pair(seq1[i - 1], seq2[j - 1]));
            int gap1_score = OPTMatrix[i - 1][j] + gap_penalty;
            int gap2_score = OPTMatrix[i][j - 1] + gap_penalty;
            OPTMatrix[i][j] = max(0, max(match_score, max(gap1_score, gap2_score)));

            // Update the maximum score and its position
            if (OPTMatrix[i][j] > max_score) {
                max_score = OPTMatrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Perform a traceback operation to locate the sequences that have been aligned
    vector<char> alignment_first_sequence, alignment_second_sequence;
    int i = max_i, j = max_j;
    while (i > 0 && j > 0 && OPTMatrix[i][j] > 0) {
        if (OPTMatrix[i][j] == OPTMatrix[i - 1][j - 1] + matrix.at(make_pair(seq1[i - 1], seq2[j - 1]))) {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), seq1[i - 1]);
            alignment_second_sequence.insert(alignment_second_sequence.begin(), seq2[j - 1]);
            i -= 1;
            j -= 1;
        } else if (OPTMatrix[i][j] == OPTMatrix[i - 1][j] + gap_penalty) {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), seq1[i - 1]);
            alignment_second_sequence.insert(alignment_second_sequence.begin(), '-');
            i -= 1;
        } else {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), '-');
            alignment_second_sequence.insert(alignment_second_sequence.begin(), seq2[j - 1]);
            j -= 1;
        }
    }

    int total_score = max_score;

    // Print the aligned sequences and total score
    cout << "Alignment of Sequence 1: ";
    for (char c : alignment_first_sequence) {
        cout << c;
    }
    cout << endl;

    cout << "Alignment of Sequence 2: ";
    for (char c : alignment_second_sequence) {
        cout << c;
    }
    cout << endl;

    // Print the OPT Matrix
    cout << "OPTMatrix: " << endl;
    for (const std::vector<int>& row : OPTMatrix) {
        for (int value : row) {
            std::cout << value << ", ";
        }
        std::cout << std::endl;
    }

    // Print optimal score
    cout << "Optimal alignment Score: " << total_score << endl;
}

// Implementation of the semi_global_alignment

void semi_global_alignment(const vector<char>& seq1, const vector<char>& seq2,
                           const map<pair<char, char>, int>& matrix, int gap_penalty) {
    int m = seq1.size();
    int n = seq2.size();

    // Initialize the table
    vector<vector<int>> OPTMatrix(m + 1, vector<int>(n + 1, 0));

    // Fill in the table
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match_score = OPTMatrix[i - 1][j - 1] + matrix.at(make_pair(seq1[i - 1], seq2[j - 1]));
            int gap1_score = OPTMatrix[i - 1][j] + gap_penalty;
            int gap2_score = OPTMatrix[i][j - 1] + gap_penalty;
            OPTMatrix[i][j] = max({0, match_score, gap1_score, gap2_score});
        }
    }

    // Find the maximum score in the last row or last column (semi-global alignment)
    int max_score = 0;
    int max_i = 0, max_j = 0;

    for (int i = 0; i <= m; ++i) {
        if (OPTMatrix[i][n] > max_score) {
            max_score = OPTMatrix[i][n];
            max_i = i;
            max_j = n;
        }
    }

    for (int j = 0; j <= n; ++j) {
        if (OPTMatrix[m][j] > max_score) {
            max_score = OPTMatrix[m][j];
            max_i = m;
            max_j = j;
        }
    }

    // Perform a traceback operation to locate the sequences that have been aligned
    vector<char> alignment_first_sequence, alignment_second_sequence;
    int i = max_i, j = max_j;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && OPTMatrix[i][j] == OPTMatrix[i - 1][j - 1] + matrix.at(make_pair(seq1[i - 1], seq2[j - 1]))) {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), seq1[i - 1]);
            alignment_second_sequence.insert(alignment_second_sequence.begin(), seq2[j - 1]);
            i -= 1;
            j -= 1;
        } else if (i > 0 && OPTMatrix[i][j] == OPTMatrix[i - 1][j] + gap_penalty) {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), seq1[i - 1]);
            alignment_second_sequence.insert(alignment_second_sequence.begin(), '-');
            i -= 1;
        } else {
            alignment_first_sequence.insert(alignment_first_sequence.begin(), '-');
            alignment_second_sequence.insert(alignment_second_sequence.begin(), seq2[j - 1]);
            j -= 1;
        }
    }

    int total_score = max_score;

    // Print the aligned sequences
    cout << "Alignment of Sequence 1: ";
    for (char c : alignment_first_sequence) {
        cout << c;
    }
    cout << endl;

    cout << "Alignment of Sequence 2: ";
    for (char c : alignment_second_sequence) {
        cout << c;
    }
    cout << endl;

    // Print the OPT Matrix
    cout << "OPTMatrix: " << endl;
    for (const std::vector<int>& row : OPTMatrix) {
        for (int value : row) {
            std::cout << value << ", ";
        }
        std::cout << std::endl;
    }

    // Print optimal score
    cout << "Optimal alignment Score: " << total_score << endl;
}

int main() {
    string sequence1, sequence2, substitution_matrix;
    char chosenalignment;
    int gap_penalty;

    /*
    Prompt the user to enter the two sequence files and substitution matrix file
    */

    cout << "Please enter amino acid sequence file 1 (with .txt): ";
    cin >> sequence1;

    cout << "Please enter amino acid sequence file 2 (with .txt): ";
    cin >> sequence2;

    cout << "Please enter substitution matrix file (with .txt): ";
    cin >> substitution_matrix;

    vector<char> testseq1 = parse_sequence_file(sequence1);
    vector<char> testseq2 = parse_sequence_file(sequence2);
    map<pair<char, char>, int> amino_Sub_matrix = parse_matrix_file(substitution_matrix);

    /*
    Prompt the user to choose the alignment algorithm
    */

    cout << "Please Select Alignment Type" << endl;
    cout << "G. Global Alignment" << endl;
    cout << "L. Local Alignment" << endl;
    cout << "S. Semi Global Alignment" << endl;

    cout << "Please Enter your choice of Alignment (G/L/S): ";
    cin >> chosenalignment;

    // Convert the entered character to uppercase if it is lowercase
    if (islower(chosenalignment)) {
        chosenalignment = toupper(chosenalignment);
    }

    switch (chosenalignment) {
        case 'G':
            cout << "You have chosen Global Alignment" << endl;
            cout << "Please enter gap penalty: ";
            cin >> gap_penalty;
            global_alignment(testseq1, testseq2, amino_Sub_matrix, gap_penalty);
            break;
        case 'L':
            cout << "You selected Local Alignment" << endl;
            cout << "Please enter gap penalty: ";
            cin >> gap_penalty;
            local_alignment(testseq1, testseq2, amino_Sub_matrix, gap_penalty);
            break;
        case 'S':
            cout << "You selected Semi Global Alignment" << endl;
            cout << "Please enter gap penalty: ";
            cin >> gap_penalty;
            semi_global_alignment(testseq1, testseq2, amino_Sub_matrix, gap_penalty);
            break;
        default:
            cout << "Invalid choice. Please enter a valid option (G/L/S)." << endl;
    }

    return 0;
}

