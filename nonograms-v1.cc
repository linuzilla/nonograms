#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cassert>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

enum CeilEnum { UNKNOWN = 0, CROSS = 1, FILLED = 2 };

class Utilities {
public:

    // 把字串拆成陣列
    static vector<int>* toIntArray(const string& s) {
    	vector<int>	*v = new vector<int>();
    	std::stringstream ss(s);
    	int input;

        while (ss >> input) {
            v->push_back(input);
        }
    
    	return v;
    }
};

class Brick {
protected:
    const int width;
    int minpos;
    int maxpos;
    bool solved;
public:
    Brick(const int w) : width(w) {
        solved = false;
        minpos = -1;
        maxpos = 999999;
    }
    
    int getWidth() const {
        return width;
    }
    
    void setMinpos(const int minpos) {
        this->minpos = this->minpos > minpos ? this->minpos : minpos; 
    }
    void setMaxpos(const int maxpos) {
        this->maxpos = this->maxpos < maxpos ? this->maxpos : maxpos;
    }
    
    int getMinpos() const {
        return minpos;
    }
    
    int getMaxpos() const {
        return maxpos;
    }
    
    bool isSatisfied() const {
        return this->maxpos - this->minpos == this->width - 1;
    }
    
    bool isSolved() const {
        return solved;
    }
    
    void markSolved(const int from) {
        solved = true;
        minpos = from;
        maxpos = from + width - 1;
    }
    
    void print() {
        cout << "{" << width << ":" << minpos << "-" << maxpos << "}";
    }
};

class LineData {
protected:
    const int row;
    const int column;
    const int numOfData;
    const int ceilSize;
    Brick       **bricks;
    CeilEnum   *ceils;
    int min;
    int firstUnsolved;
    int lastUnsolved;
    int leftMargin;
    int rightMargin;
    bool solved;
public:
    LineData(const int row, const int column, const int numOfData, const int ceilSize, vector<int> *data)
            : row(row), column(column), numOfData(numOfData), ceilSize(ceilSize) {
        solved = false;
        ceils = new CeilEnum[ceilSize];
        leftMargin = 0;
        rightMargin = ceilSize;
        
        if (numOfData == 0) {
            min = 0;
            for (int i = 0; i < ceilSize; i++) {
                ceils[i] = CROSS;
            }
            solved = true;
        } else {
            bricks = new Brick*[numOfData];
            
            firstUnsolved = 0;
            lastUnsolved = numOfData - 1;
            
            min = numOfData - 1;
            
            int minpos = 0;
            
            for (int i = 0; i < numOfData; i++) {
                bricks[i] = new Brick(data->at(i));
                bricks[i]->setMinpos(minpos);
                minpos += bricks[i]->getWidth() + 1;
                min += data->at(i);
            }
            
            int maxpos = ceilSize - 1;
            
            for (int i = numOfData - 1; i >= 0; i--) {
                bricks[i]->setMaxpos(maxpos);
                maxpos -= bricks[i]->getWidth() + 1;
            }
            
            if (min == ceilSize) {
                solved = true;
                
                int j = 0;
                int k = 0;
                
                for (int i = 0; i < ceilSize; i++) {
                    if (k++ < data->at(j)) {
                        ceils[i] = FILLED;
                    } else {
                        j++;
                        k = 0;
                        ceils[i] = CROSS;
                    }
                }
            } else {
                for (int i = 0; i < ceilSize; i++) {
                    ceils[i] = UNKNOWN;
                }
            }
        }
        
        data->clear();
        delete data;
    };
    
    void print() const {
        cerr << "{";
        for (int i = 0; i < numOfData; i++) {
            cerr << " " <<  bricks[i]->getWidth();
            if (bricks[i]->isSolved()) cerr << "s";
        }
        cerr << " } [" << firstUnsolved << "," << lastUnsolved << "]" << (solved ? " solved" : "");
    }
    
    void display() const {
        cerr << "{ ";
        for (int i = 0; i < ceilSize; i++) {
            if (i % 10 == 0) cerr << "|";
            switch (ceils[i]) {
            case UNKNOWN:
                cerr << " ";
                break;
            case CROSS:
                cerr << ".";
                break;
            case FILLED:
                cerr << "#";
                break;
            }
        }
        cerr << " }";
    }
    
    int getMin() const {
        return min;
    }
    
    bool isSolved() const {
        return solved;
    }
    
    void markSolved() {
        solved = true;
        
        for (int i = 0; i < ceilSize; i++) {
            if (ceils[i] == UNKNOWN)
                ceils[i] = CROSS;
        }
    }
    
    bool checkIsSolved() {
        if (! solved) {
            for (int i = 0; i < ceilSize; i++)
                if (ceils[i] == UNKNOWN)
                    return false;
            solved = true;
        }
        return solved;
    }
    
    CeilEnum * getCeils() {
        return ceils;
    }
    
    void makeItSolved(const int i, const int from) {
        cerr << "MakeItSolved " << i << " from:" << from << "   ";
        this->print(); cerr << endl;
        this->display(); cerr << endl;
        for (int j = 0, k = from; j < bricks[i]->getWidth(); k++, j++) {
            switch (ceils[k]) {
            case UNKNOWN:
                ceils[k] = FILLED;
                break;
            case FILLED:
                break;
            case CROSS:
                assert(false);
                break;
            }
        }

        bricks[i]->markSolved(from);
        
        if (from > 0) {
            assert(ceils[from - 1] != FILLED);
            ceils[from - 1] = CROSS;
        }
        
        int r = from + bricks[i]->getWidth();
        
        if (r < ceilSize - 1) {
            // cout << "from=" << from << ", width=" << bricks[i]->getWidth() << endl;
            if (ceils[r] == FILLED) {
                cerr << "r=" << r << endl;
            }
            assert(ceils[r] != FILLED);
            ceils[r] = CROSS;
        }
        
        if (i > 0 && ! bricks[i - 1]->isSolved()) {
            bricks[i - 1]->setMaxpos(from - 2);
        }
        if (i < numOfData - 1 && ! bricks[i + 1]->isSolved()) {
            bricks[i + 1]->setMinpos(r + 1);
        }
    }
    
    void findBySegment() {
        if (solved) return;
        
        int numOfSegments = 0;
        bool infilled = false;
        
        for (int i = 0; i < ceilSize; i++) {
            if (infilled) {
                if (ceils[i] == CROSS) {
                    infilled = false;
                }
            } else if (ceils[i] == FILLED) {
                infilled = true;
                numOfSegments++;
            }
        }
        
        if (numOfSegments == numOfData) {
            cerr << "(" << row << "," << column << ") seg= " << numOfSegments << endl;
            
            int seg = 0;
            infilled = false;
            int crosspos = 0;
            int firstFilled = 0;
            int lastFilled = -1;
            
            for (int i = 0; i < ceilSize; i++) {
                if (infilled) {
                    switch(ceils[i]) {
                    case CROSS: {
                        infilled = false;
                        
                        for (int j = firstFilled + 1; j < lastFilled; j++) {
                            if (ceils[j] == UNKNOWN) ceils[j] = FILLED;
                        }
                        
                        int width = bricks[seg]->getWidth();
                        bricks[seg]->setMinpos(lastFilled - width + 1);
                        bricks[seg]->setMaxpos(firstFilled + width - 1);

                        int minpos = bricks[seg]->getMinpos();
                        int maxpos = bricks[seg]->getMaxpos();
                    
                        for (int j = crosspos; j < minpos; j++) {
                            if (ceils[j] == FILLED) {
                                this->print(); cerr << endl;
                                this->display(); cerr << endl;
                                cerr << "(min=" << bricks[seg]->getMinpos() << ",max=" << bricks[seg]->getMaxpos() << ",width="
                                    << bricks[seg]->getWidth() << ")" << endl;
                                cerr << "seg " << seg << " fill cross (" << row << "," << column << ") from " << crosspos << " to " << minpos << ", fill " << j << endl;
                            }
                            assert (ceils[j] != FILLED);
                            if (ceils[j] == UNKNOWN) {
                                // cerr << "@@@ (" << row << "," << column << ") from " << crosspos << " to " << minpos << ", fill " << j << endl;
                                ceils[j] = CROSS;
                            }
                        }

                        // this->display(); cerr << endl;
                        
                        
                        crosspos = i + 1;

                        if (firstFilled > 0 && ceils[firstFilled - 1] == CROSS) {
                            // cerr << "[+] Solve (" << row << "," << column << ") , seg/nsg = " << seg <<  "/" << numOfData << ", first " << firstFilled << ", width " << bricks[seg]->getWidth() << endl;
                            makeItSolved(seg, firstFilled);
                        } else if (lastFilled == i - 1) {
                            makeItSolved(seg, i - bricks[seg]->getWidth());
                        }
                        
                        if (lastFilled - firstFilled + 1 == bricks[seg]->getWidth()) {
                            cerr << "(" << row << "," << column << ") "<< "last " << lastFilled << ", first " << firstFilled << ", vs " << bricks[seg]->getWidth() << endl;
                            if (firstFilled > 0) {
                                assert (ceils[firstFilled - 1] != FILLED);
                                ceils[firstFilled - 1] = CROSS;
                            }
                            if (lastFilled < ceilSize - 1) {
                                assert (ceils[lastFilled + 1] != FILLED);
                                ceils[lastFilled + 1] = CROSS;
                            }
                        }
                        
                        seg++;
                        break;
                    }

                    case FILLED:
                        lastFilled = i;
                        break;
                        
                    case UNKNOWN:
                        break;
                    }
                } else if (ceils[i] == FILLED) {
                    infilled = true;
                    firstFilled = i;
                }
            }
        }
    }
    
    bool checkComplete() {
        if (solved) return true;
        for (int i = 0; i < numOfData; i++) {
            if (! bricks[i]->isSolved()) {
                if (bricks[i]->isSatisfied()) {
                    this->makeItSolved(i, bricks[i]->getMinpos());
                } else {
                    return false;
                }
            }
        }
        this->markSolved();
        return true;
    }
    void findPossible() {
        if (solved) return;
        
        bool advance_firstUnsolved = false;
        bool advance_lastUnsolved = false;
        
        for (int i = firstUnsolved; i <= lastUnsolved; i++) {
            if (bricks[i]->isSolved()) continue;
            
            int width = bricks[i]->getWidth();
            
            bricks[i]->setMinpos(leftMargin);
            bricks[i]->setMaxpos(rightMargin);
            
            int minpos = bricks[i]->getMinpos();
            int maxpos = bricks[i]->getMaxpos();
            
            //minpos = minpos > leftMargin ? minpos : leftMargin;
            // maxpos = maxpos < rightMargin ? maxpos : rightMargin;
            
            if (i > 0 && bricks[i - 1]->isSolved()) {
                int lm = bricks[i - 1]->getMaxpos() + 1;
                minpos = minpos > lm ? minpos : lm;
                // cerr << "brick " << ( i - 1 ) << " max=" << bricks[i-1]->getMaxpos() << ", minpos=" << minpos << endl;
            }
            
            for (int j = minpos; j < maxpos; j++) {
                if (ceils[j] == CROSS) {
                    minpos = j + 1;
                } else if (i == firstUnsolved && ceils[j] == FILLED) {
                    cerr << "Try Make it solve from " << __LINE__ << " i=" << i << ", from=" << j << "(width=" << width << ":" << minpos << "," << maxpos << ")" << endl;
                    makeItSolved(i, j);
                    advance_firstUnsolved = true;
                    leftMargin = i + width + 1;
                    break;
                } else {
                    break;
                }
            }
            
            for (int j = maxpos; j > minpos; j--) {
                if (ceils[j] == CROSS) {
                    maxpos = j - 1;
                } else if (i == lastUnsolved && ceils[j] == FILLED) {
                    cerr << "Try Make it solve from " << __LINE__ << " i=" << i << ", from=" << (j - width + 1) << "(width=" << width << ":" << minpos << "," << maxpos << ")" << endl;
                    makeItSolved(i, j - width + 1);
                    rightMargin = j - width;
                    advance_lastUnsolved = true;
                    break;
                } else {
                    break;
                }
            }
            
            bricks[i]->setMinpos(minpos);
            bricks[i]->setMaxpos(maxpos);
            
            if (bricks[i]->isSolved()) {
                continue;
            } else if (width == maxpos - minpos + 1) {
                cerr << "Try Make it solve from " << __LINE__ << " i=" << i << ", from=" << minpos << "(width=" << width << ":min=" << minpos << ",max=" << maxpos << ")" << endl;
                makeItSolved(i, minpos);
                
                if (i == firstUnsolved) {
                    for (int j = rightMargin; j < minpos; j++) {
                        assert (ceils[j] != FILLED);
                        // if (ceils[j] == UNKNOWN)
                        ceils[j] = CROSS;
                    }
                    advance_firstUnsolved = true;
                } else if (i == lastUnsolved) {
                    cerr << "Line " << __LINE__ << " i=" << i << ", from=" << minpos << "(" << width << ":" << minpos << "," << maxpos << ")" << endl;
                    for (int j = maxpos + 1; j < leftMargin; j++) {
                        assert (ceils[j] != FILLED);
                        ceils[j] = CROSS;
                    }
                    advance_lastUnsolved = true;
                }
            } else if (width > maxpos - minpos + 1) {
                cout << "width=" << width << ", min=" << minpos << ", max=" << maxpos << endl;
                assert(false);
            } else if (i == firstUnsolved && i == lastUnsolved) {
                bool infilled = false;
                int count = 0;
                int firstpos = -1;
                
                for (int j = minpos; j <= maxpos; j++) {
                    if (infilled) {
                        if (ceils[j] == FILLED) {
                            count++;
                        } else {
                            break;
                        }
                    } else if (ceils[j] == FILLED) {
                        infilled = true;
                        firstpos = j;
                        count = 1;
                    }
                }
                
                if (infilled && count == width) {
                    // cout << "(" << row << "," << column << ") find " << i << " at " << firstpos << ", width=" << width << endl;
                    cerr << "Try Make it solve from " << __LINE__ << endl;
                    makeItSolved(i, firstpos);
                    advance_firstUnsolved = true;
                }
            }
            
            if (maxpos -  minpos + 1 < width + width) {
                for (int j = maxpos - width + 1; j < minpos + width; j++) {
                    if (ceils[j] == UNKNOWN) {
                        ceils[j] = FILLED;
                    } else if (ceils[j] == CROSS) {
                        assert (false);
                    }
                }
            }
        }
        
        if (firstUnsolved > lastUnsolved) {
            cout << "solved (" << row << "," << column << "): ";
            this->print();
            for (int i = 0; i < ceilSize; i++) {
                if (ceils[i] == UNKNOWN) {
                    cout << " " << i;
                    ceils[i] = CROSS;
                }
            }
            this->markSolved();
            cout << endl;
        }
        
        if (advance_firstUnsolved) firstUnsolved++;
        if (advance_lastUnsolved) lastUnsolved--;
    }
};

class Nonogram {
protected:
    const int row;
    const int col;
    LineData **rowdata;
    LineData **coldata;
    CeilEnum **ceils;
public:
    Nonogram(const int row, const int col) : row(row), col(col) {
        rowdata = new LineData*[row];
        coldata = new LineData*[col];
        ceils = new CeilEnum*[row];
        for (int i = 0; i < row; i++)
            ceils[i] = new CeilEnum[col];
            
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                ceils[i][j] = UNKNOWN;
    }
    
    void setRow(const int r, vector<int> *data) {
        rowdata[r] = new LineData(r, -1, data->size(), col, data);
    }
    
    void setCol(const int c, vector<int> *data) {
        coldata[c] = new LineData(-1, c, data->size(), row, data);
    }
    
    void print() const {
        for (int i = 0; i < row; i++) {
            cout << std::setw(2) << i << ": ";
            for (int j = 0; j < col; j++) {
                switch (ceils[i][j]) {
                case UNKNOWN:
                    cout << " ";
                    break;
                case FILLED:
                    cout << "#";
                    break;
                case CROSS:
                    cout << ".";
                    break;
                }
            }
            cout << "  : ";
            rowdata[i]->print();
            cout << endl;
        }
    }
    
    bool findAndSync(const int r, const int c) {
        LineData *linedata;
        int max;
        int x = 0;
        int y = 0;
        int dx = 0;
        int dy = 0;
        
        if (r != -1) {
            y = r;
            dx = 1;
            max = col;
            linedata = rowdata[r];
        } else if (c != -1) {
            x = c;
            dy = 1;
            max = row;
            linedata = coldata[c];
        }
        
        linedata->findPossible();
        linedata->findBySegment();
        linedata->checkComplete();
        
        CeilEnum *ceilsData = linedata->getCeils();
        
        for (int i = 0; i < max; i++, x += dx, y += dy) {
            if (ceilsData[i] != ceils[y][x]) {
                if (ceilsData[i] == UNKNOWN) {
                    ceilsData[i] = ceils[y][x]; 
                } else if (ceils[y][x] == UNKNOWN) {
                    ceils[y][x] = ceilsData[i];
                } else {
                    cout << "i=" << i << ",x=" << x << ",y=" << y << ", dx=" << dx << ", dy=" << dy << endl;
                    cout << "ceils[" << y << "][" << x << "]=" << ceilsData[i] << ", " << ceils[y][x] << endl;
                    assert(false);
                }
            }
            // this->print();
        }
        
        return linedata->checkIsSolved();
    }
    
    int trySolve(const int r, const int c) {
        return findAndSync(r, c) ? 1 : 0;
    }
    
    void solve() {
        for (int i = 1; i < 16; i++) {
            int count = 0;
            for (int r = 0; r < row; r++)
                count += trySolve(r, -1);
    
            for (int c = 0; c < col; c++)
                count += trySolve(-1, c);
                
            cout << "Solved Count=" << count << endl;
        }
            
        // cout << "Solved Count=" << count << endl;
    }
};

Nonogram* inputReader (const char *filename) {
    Nonogram *nonogram = NULL;
    std::ifstream ifs(filename);
    string line;
    
    if (ifs.is_open()) {
        if (std::getline (ifs, line)) {
            vector<int> *list = Utilities::toIntArray(line);
            int row = list->at(0);
            int col = list->at(1);
            list->clear();
            delete list;
            
            nonogram = new Nonogram(row, col);
            
            for (int r = 0; r < row; r++) {
                std::getline(ifs, line);
                nonogram->setRow(r, Utilities::toIntArray(line));
            }
            
            for (int c = 0; c < col; c++) {
                std::getline(ifs, line);
                nonogram->setCol(c, Utilities::toIntArray(line));
            }
        }
        ifs.close();
    }
    return nonogram;
}

int main (int argc, char *argv[]) {
    if (argc > 1) {
        Nonogram *nonogram = inputReader (argv[1]);
        nonogram->solve();
        nonogram->print();
    } else {
        cout << "usage: nonograms inputfile" << endl;
    }
    return 0;
}