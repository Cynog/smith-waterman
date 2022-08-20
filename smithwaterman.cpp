#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

int smithwaterman(string seq1, string seq2, int match, int mismatch, int gap);

int main(int argc, char *argv[]) {
    fstream file;
    file.open("seq1.txt", ios::in);
    string seq1;
    getline(file, seq1);
    file.close();
    file.open("seq2.txt", ios::in);
    string seq2;
    getline(file, seq2);
    file.close();

    int score = smithwaterman(seq1, seq2, 1, -1, -1);

    cout << "score = " << score << endl;

    // return
    return 0;
}

int smithwaterman(string seq1, string seq2, int match, int mismatch, int gap) {
    int len1 = seq1.length();
    int len2 = seq2.length();

    vector<vector<int>> scores(len1 + 1, vector<int>(len2 + 1, 0));
    vector<vector<char>> paths(len1 + 1, vector<char>(len2 + 1, 'o'));

    // fill smith waterman table
    for (int k = 2; k <= len1 + len2; k++) {
#pragma omp parallel for
        for (int i = 1; i <= len1; i++) {
            int j = k - i;
            if (j < 1 or j > len2) continue;

            int diag;
            if (seq1[i - 1] == seq2[j - 1])
                diag = scores[i - 1][j - 1] + match;
            else
                diag = scores[i - 1][j - 1] + mismatch;

            int down = scores[i - 1][j] + gap;
            int right = scores[i][j - 1] + gap;

            int max = 0;
            if (diag > max) {
                max = diag;
                paths[i][j] = 'd';
            }
            if (down > max) {
                max = down;
                paths[i][j] = 'u';
            }
            if (right > max) {
                max = right;
                paths[i][j] = 'l';
            }

            scores[i][j] = max;
        }
    }

    // search for highest score
    int i_max = 0;
    int j_max = 0;
    int score = 0;
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            // cout << scores[i][j] << " ";
            if (scores[i][j] > score) {
                score = scores[i][j];
                i_max = i;
                j_max = j;
            }
        }
        // cout << endl;
    }

    string seq1a = "";
    string seq2a = "";
    int i = i_max;
    int j = j_max;
    while (scores[i][j] != 0) {
        if (paths[i][j] == 'd') {
            seq1a.insert(0, 1, seq1[i - 1]);
            seq2a.insert(0, 1, seq2[j - 1]);
            i -= 1;
            j -= 1;
        } else if (paths[i][j] == 'u') {
            seq1a.insert(0, 1, seq1[i - 1]);
            seq2a.insert(0, 1, '-');
            i -= 1;
        } else if (paths[i][j] == 'l') {
            seq1a.insert(0, 1, '-');
            seq2a.insert(0, 1, seq2[j - 1]);
            j -= 1;
        }
    }

    cout << seq1a << endl;
    cout << seq2a << endl;

    // return local alignment score
    return score;
}