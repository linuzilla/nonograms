#include <iostream>
#include <sstream>
#include <vector>
#include <stack>
#include <cstdlib>
#include <iomanip>
#include <cassert>

static bool trace = false;

#if ENABLE_CTRL_C_TRACE
#include <signal.h>
void ctrlc_handler(int signo) {
    trace = ! trace;
}
#endif


using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


enum CeilEnum { UNKNOWN = 0, CROSS = 1, FILLED = 2 };

class Utilities {
public:
    static vector<int>* toIntArray(const string& s) { // 把字串拆成陣列
        vector<int>    *v = new vector<int>();
        std::stringstream ss(s); // 把字串改成 stream
        int input;

        while (ss >> input) // 從 string stream 取資料
            v->push_back(input); // 塞到 vector 裡
    
        return v;
    }
};

class Bead {
protected:
    const int width;
    int minpos;
    int maxpos;
    bool solved;
public:
    Bead(const int w) : width(w) {
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
};

class BeadsState {
protected:
    const int numOfBeads;
    int maxValue;
    int start;
    int last;
    int maxAllowed;
    int stateValue;
    int *states;
    bool initialized;

public:
    BeadsState(const int num, const int max) : numOfBeads(num), maxValue(max) {
        states = new int[num];

        for (int i = 0; i < numOfBeads; i++) {
            states[i] = 0;
        }

    start = 0;
        stateValue = 0;
        
        maxAllowed = maxValue;
        initialized = false;
        last = numOfBeads;
    };

    BeadsState(BeadsState& from) : numOfBeads(from.numOfBeads), maxValue(from.maxValue) {
        states = new int[numOfBeads];

        for (int i = 0; i < numOfBeads; i++) {
            states[i] = from.states[i];
        }
        stateValue = from.stateValue;
        maxAllowed = from.maxAllowed;
        initialized = from.initialized;
        last = from.last;
    start = from.start;
    }

    ~BeadsState() {
        delete states;
    }

    BeadsState* duplicate() {
        return new BeadsState(*this);
    }

    int getBead(const int i) const {
        return states[i];
    }
    
    void setMaxValue (const int v) {
        maxValue = maxValue < v ? maxValue : v;
    }
    
    void setLast (const int lastv) {
        last = lastv;
    }
    
    int getLast () const {
        return last;
    }

    void setStart(const int start) {
        this->start = start;
    }

    void setBead(const int i, const int v) {
        states[i] = v;
        initialized = false;
    }
    
    bool nextState() {    
        if (initialized == false) {
            initialized = true;
            maxAllowed = maxValue;

            for (int i = 0; i < last; i++) {
                maxAllowed -= states[i];
            }
            return true;
        } else {
            for (int i = last - 1; i >= start; i--) {
                if (maxAllowed > 0) {
                    states[i]++;
                    maxAllowed--;

                    return true;
                } else {
                    maxAllowed += states[i];
                    states[i] = 0;
                }
            }
        }
        
        return false;
    }

    friend std::ostream& operator<<(std::ostream& os, BeadsState& bs);
};

std::ostream& operator<<(std::ostream& os, BeadsState& bs) {
    os << "[";
    for (int i = 0; i < bs.numOfBeads; i++) {
        os << " " << bs.states[i];
    }
    os << " ]";
    return os;
}

class Axis {
protected:
    const int numOfBeads;
    const int row;
    const int column;
    const int ceilSize;
    Bead           **beads;
    CeilEnum       *ceils;
    CeilEnum    *initialguess_ceils;
    BeadsState    *state;
    std::stack<BeadsState *> stateStack;
    int min;
    int firstUnsolved;
    int lastUnsolved;
    int leftMargin;
    int rightMargin;
    bool solved;
public:
    Axis(const int row, const int column, const int numOfBeads, const int ceilSize, vector<int> *data)
            : numOfBeads(numOfBeads), row(row), column(column), ceilSize(ceilSize) {
        solved = false;
        ceils = new CeilEnum[ceilSize];
        initialguess_ceils = new CeilEnum[ceilSize];
        
        leftMargin = 0;
        rightMargin = ceilSize;
        
        if (numOfBeads == 0) {
            beads = NULL;
            state = NULL;
            min = 0;
            for (int i = 0; i < ceilSize; i++) {
                ceils[i] = CROSS;
            }
            solved = true;
        } else {
            beads = new Bead*[numOfBeads];
            
            int maxValue = ceilSize - numOfBeads + 1;
            firstUnsolved = 0;
            lastUnsolved = numOfBeads - 1;
            
            min = numOfBeads - 1;
            
            int minpos = 0;
        
            for (int i = 0; i < numOfBeads; i++) {
                beads[i] = new Bead(data->at(i));
                beads[i]->setMinpos(minpos);
                minpos += beads[i]->getWidth() + 1;
                min += data->at(i);
                maxValue -= beads[i]->getWidth();
            }
            state = new BeadsState(numOfBeads, maxValue);
            
            int maxpos = ceilSize - 1;
            
            for (int i = numOfBeads - 1; i >= 0; i--) {
                beads[i]->setMaxpos(maxpos);
                maxpos -= beads[i]->getWidth() + 1;
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
    }
    
    ~Axis() {
        if (numOfBeads > 0) {
        }
    }

    BeadsState *getState() {
        return state;
    }

    int getNumOfBeads() const {
        return numOfBeads;
    }

    bool isSolved() const {
        return solved;
    }
    
    bool nextState() {
        // BeadsState *newState = state->duplicate();
        // cerr << "Row=" << row << " ";
        bool ns = state->nextState();
        // if (trace)
        //     cerr << "Row:" << row << ", State:" << *state << endl;
        return ns;
    }
    
    void pushState() {
        stateStack.push(state->duplicate());
    }

    void popState() {
        delete state;
        state = stateStack.top();
        stateStack.pop();
    }
    
    CeilEnum * getCeils() {
        return ceils;
    }
    
    void markSolved() {
        solved = true;
        
        for (int i = 0; i < ceilSize; i++) {
            if (ceils[i] == UNKNOWN)
                ceils[i] = CROSS;
        }
    }

    bool updateCeilsByState() {
        for (int i = 0; i < ceilSize; i++) {
            ceils[i] = CROSS;
        }

        int pos = 0;
        int last = state->getLast();

        for (int i = 0; i < last; i++) {
            assert (state->getBead(i) >= 0);
            pos += state->getBead(i);
            for (int j = 0; j < beads[i]->getWidth(); j++) {
                ceils[pos++] = FILLED;
            }

            pos++;
        }
        
        for (int i = 0; i < pos; i++)
            if (initialguess_ceils[i] != UNKNOWN && initialguess_ceils[i] != ceils[i])
                return false;
        
        for (int i = pos; i < ceilSize; i++)
            ceils[i] = initialguess_ceils[i] == UNKNOWN ? CROSS : initialguess_ceils[i];

        return true;
    }
    
    void print(std::ostream& os) const {
        os << "{";
        for (int i = 0; i < numOfBeads; i++) {
            os << " " <<  beads[i]->getWidth();
            if (beads[i]->isSolved()) os << "s";
        }
        os << " } [" << firstUnsolved << "," << lastUnsolved << "]" << (solved ? " solved" : "");
    }
    
    void display(std::ostream& os) const {
        for (int i = 0; i < ceilSize; i++) {
            // if (trace && (i % 10 == 0)) os << "|";
            if (ceils[i] != initialguess_ceils[i] && trace) {
                switch (initialguess_ceils[i]) {
                case UNKNOWN:
                    if (ceils[i] == FILLED) {
                        os << "@";
                    } else {
                        os << ",";
                    }
                    break;
                case CROSS:
                    os << ".";
                    break;
                case FILLED:
                    os << "#";
                    break;
                }
            } else {
                switch (ceils[i]) {
                case UNKNOWN:
                    os << " ";
                    break;
                case CROSS:
                    os << ".";
                    break;
                case FILLED:
                    os << "#";
                    break;
                }
            }
        }
        // cerr << " }";
    }
    
    bool checkComplete() {
        if (solved) return true;
        for (int i = 0; i < numOfBeads; i++) {
            if (! beads[i]->isSolved()) {
                if (beads[i]->isSatisfied()) {
                    this->makeItSolved(i, beads[i]->getMinpos());
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
            if (beads[i]->isSolved()) continue;
            
            int width = beads[i]->getWidth();
            
            beads[i]->setMinpos(leftMargin);
            beads[i]->setMaxpos(rightMargin);
            
            int minpos = beads[i]->getMinpos();
            int maxpos = beads[i]->getMaxpos();
            
            //minpos = minpos > leftMargin ? minpos : leftMargin;
            // maxpos = maxpos < rightMargin ? maxpos : rightMargin;
            
            if (i > 0 && beads[i - 1]->isSolved()) {
                int lm = beads[i - 1]->getMaxpos() + 1;
                minpos = minpos > lm ? minpos : lm;
            }
            
            for (int j = minpos; j < maxpos; j++) {
                if (ceils[j] == CROSS) {
                    minpos = j + 1;
                } else if (i == firstUnsolved && ceils[j] == FILLED) {
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
                    makeItSolved(i, j - width + 1);
                    rightMargin = j - width;
                    advance_lastUnsolved = true;
                    break;
                } else {
                    break;
                }
            }
            
            beads[i]->setMinpos(minpos);
            beads[i]->setMaxpos(maxpos);
            
            if (beads[i]->isSolved()) {
                continue;
            } else if (width == maxpos - minpos + 1) {
                makeItSolved(i, minpos);
                
                if (i == firstUnsolved) {
                    for (int j = rightMargin; j < minpos; j++) {
                        assert (ceils[j] != FILLED);
                        ceils[j] = CROSS;
                    }
                    advance_firstUnsolved = true;
                } else if (i == lastUnsolved) {
                    for (int j = maxpos + 1; j < leftMargin; j++) {
                        assert (ceils[j] != FILLED);
                        ceils[j] = CROSS;
                    }
                    advance_lastUnsolved = true;
                }
            } else if (width > maxpos - minpos + 1) {
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
            for (int i = 0; i < ceilSize; i++) {
                if (ceils[i] == UNKNOWN) {
                    ceils[i] = CROSS;
                }
            }
            this->markSolved();
        }
        
        if (advance_firstUnsolved) firstUnsolved++;
        if (advance_lastUnsolved) lastUnsolved--;
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
        
        if (numOfSegments == numOfBeads) {
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
                        
                        int width = beads[seg]->getWidth();
                        beads[seg]->setMinpos(lastFilled - width + 1);
                        beads[seg]->setMaxpos(firstFilled + width - 1);

                        int minpos = beads[seg]->getMinpos();
                        // int maxpos = beads[seg]->getMaxpos();
                    
                        for (int j = crosspos; j < minpos; j++) {
                            assert (ceils[j] != FILLED);
                            if (ceils[j] == UNKNOWN)
                                ceils[j] = CROSS;
                        }
                        
                        
                        crosspos = i + 1;

                        if (firstFilled > 0 && ceils[firstFilled - 1] == CROSS) {
                            makeItSolved(seg, firstFilled);
                        } else if (lastFilled == i - 1) {
                            makeItSolved(seg, i - beads[seg]->getWidth());
                        }
                        
                        if (lastFilled - firstFilled + 1 == beads[seg]->getWidth()) {
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
    
    void findByParsing() {
        if (solved) return;
        
        bool inBead = false;
        int minBeadId = 0;
        int beadLength = 0;

        for (int i = 0; i < ceilSize; i++) {
            if (inBead) {
                if (ceils[i] == FILLED) {
                    beadLength++;
                } else {
                    inBead = false;
                    // Check beed belonging
                    for (int j = minBeadId; j < numOfBeads; j++) {
                        if (beadLength > beads[j]->getWidth()) {
                            minBeadId = j + 1;
                        } else {
                            if (beads[j]->getMaxpos() < i - 1) {
                                minBeadId = j + 1;
                            } else {
                                break;
                            }
                        }
                    }
                    
                    assert (minBeadId < numOfBeads);
                    
                    if (minBeadId == numOfBeads - 1)
                        if (beadLength == beads[minBeadId]->getWidth())
                            makeItSolved(minBeadId, i - beadLength);
                }
            } else {
                if (ceils[i] == FILLED) {
                    inBead = true;
                    beadLength = 1;
                    while (beads[minBeadId]->getMaxpos() < i)
                        minBeadId++;

                    assert (beads[minBeadId]->getMinpos() <= i);
                }
            }
        }
    }
    
    void findImpossible () {
        if (solved) return;
        
        bool inBead = false;
        int minBeadId = 0;
        int beadLength = 0;
        int startPos = 0;

        for (int i = 0; i < ceilSize; i++) {
            if (inBead) {
                if (ceils[i] == FILLED || ceils[i] == UNKNOWN) {
                    beadLength++;
                } else {
                    if (beads[minBeadId]->getWidth() > beadLength) {
                        for (int j = i - 1; j >= startPos; j--) {
                            assert (ceils[i] != FILLED);
                            ceils[j] = CROSS;
                        }
                    }
                    break;
                }
            } else {
                if (ceils[i] == FILLED || ceils[i] == UNKNOWN) {
                    inBead = true;
                    beadLength = 1;
                    startPos = i;
                }
            }
        }
    }
    
    void makeItSolved(const int i, const int from) {
        for (int j = 0, k = from; j < beads[i]->getWidth(); k++, j++) {
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

        beads[i]->markSolved(from);
        
        if (from > 0) {
            assert(ceils[from - 1] != FILLED);
            ceils[from - 1] = CROSS;
        }
        
        int r = from + beads[i]->getWidth();
        
        if (r < ceilSize) {
            assert(ceils[r] != FILLED);
            ceils[r] = CROSS;
        }
        
        if (i > 0 && ! beads[i - 1]->isSolved()) {
            beads[i - 1]->setMaxpos(from - 2);
        }
        if (i < numOfBeads - 1 && ! beads[i + 1]->isSolved()) {
            beads[i + 1]->setMinpos(r + 1);
        }
    }
    
    void updateState() {
        int pos = 0;
        bool inSolved = true;
        int start = 0;
        
        for (int i = 0; i < numOfBeads; i++) {
            beads[i]->setMinpos(pos);
            state->setBead(i, beads[i]->getMinpos() - pos);
            assert (state->getBead(i) >= 0);
            pos += state->getBead(i) + beads[i]->getWidth() + 1;

            if (inSolved && beads[i]->isSolved()) {
                start = i + 1;
            } else {
                inSolved = false;
            }
        }
    
        state->setStart(start);
        
        for (int i = numOfBeads - 1; i > 0; i--) {
            if (beads[i]->isSolved()) {
                state->setLast(i);
                
                int maxValue = beads[i]->getMinpos() - i;
        
                for (int j = 0; j < i; j++) {
                    maxValue -= beads[j]->getWidth();
                }
                state->setMaxValue(maxValue);
                assert (maxValue >= 0);
            } else {
                break;
            }
        }
        
        for (int i = 0; i < ceilSize; i++) {
            initialguess_ceils[i] = ceils[i];
        }
    }
    
    bool combination_check (const int r, CeilEnum comb) {
        return CeilEnum(ceils[r] & comb) != ceils[r];
    }
    
    CeilEnum calculate_allowed_combination (const int r) {
        bool inBead = false;
        int beadId = 0;
        int beadLength = 0;

        for (int i = 0; i < r; i++) {
            if (inBead) {
                if (ceils[i] == FILLED) {
                    beadLength++;
                } else {
                    inBead = false;
                    beadId++;
                }
            } else {
                if (ceils[i] == FILLED) {
                    inBead = true;
                    beadLength = 1;
                }
            }
        }
        
        if (inBead) {
            return beadLength < beads[beadId]->getWidth() ? FILLED : CROSS; 
        } else {
            return (beadId == numOfBeads) ? CROSS : CeilEnum (CROSS | FILLED);
        }
    }

    bool checkByCeils() {
        bool inBead = false;
        int beadId = 0;
        int beadLength = 0;

        for (int i = 0; i < ceilSize; i++) {
            if (inBead) {
                if (ceils[i] == FILLED) {
                    beadLength++;
                } else {
                    inBead = false;
                    if (beads[beadId]->getWidth() != beadLength) return false;
                    beadId++;
                }
            } else {
                if (ceils[i] == FILLED) {
                    inBead = true;
                    beadLength = 1;
                    if (beadId >= numOfBeads) return false;
                }
            }
        }
        if (inBead) {
            if (beads[beadId]->getWidth() != beadLength) return false;
        }
        return true;
    }
};

class Nonogram {
protected:
    const int row;
    const int col;
    Axis **rowdata;
    Axis **coldata;
    CeilEnum **ceils;

    bool rewrite_row(const int r) {
        if (rowdata[r]->updateCeilsByState()) {
            CeilEnum *ceilsData = rowdata[r]->getCeils();

            for (int i = 0; i < col; i++)
                coldata[i]->getCeils()[r] = ceils[r][i] = ceilsData[i];

            return true;
        } else {
            return false;
        }
    }

    bool early_checker(const int r, CeilEnum *comb) {
        for (int i = 0; i < col; i++)
            if (coldata[i]->combination_check (r, comb[i]))
                return false;

        return true;
    }
    
    CeilEnum *  allowed_combination (const int r) {
        CeilEnum *pattern = new CeilEnum[col];
        for (int i = 0; i < col; i++)
            pattern[i] = coldata[i]->calculate_allowed_combination(r);
        return pattern;
    }
    
    bool findAndSync(const int r, const int c) {
        Axis *axis;
        int max ;
        int x = 0;
        int y = 0;
        int dx = 0;
        int dy = 0;
        int progress = 0;
        
        if (r != -1) {
            y = r;
            dx = 1;
            max = col;
            axis = rowdata[r];
        } else {
            x = c;
            dy = 1;
            max = row;
            axis = coldata[c];
        }
        
        axis->findPossible();
        axis->findBySegment();
        axis->findByParsing();
        axis->findImpossible();
        axis->checkComplete();
        
        CeilEnum *ceilsData = axis->getCeils();
        
        for (int i = 0; i < max; i++, x += dx, y += dy) {
            if (ceilsData[i] != ceils[y][x]) {
                if (ceilsData[i] == UNKNOWN) {
                    ceilsData[i] = ceils[y][x];
                    progress++;
                } else if (ceils[y][x] == UNKNOWN) {
                    ceils[y][x] = ceilsData[i];
                    progress++;
                } else {
                    assert(false);
                }
            }
        }

        return progress;
    }
    
    int guess_iteration() {
        int progress = 0;
        for (int r = 0; r < row; r++)
            progress += findAndSync(r, -1);
        
        for (int c = 0; c < col; c++)
            progress += findAndSync(-1, c);
        return progress;
    }
    
    void initial_guess() {
        while (guess_iteration() > 0)
            ;
        
        for (int r = 0; r < row; r++)
            rowdata[r]->updateState();
    }
    
public:
    static Nonogram* inputReader ();
    
    Nonogram(const int row, const int col) : row(row), col(col) {
        rowdata = new Axis*[row];
        coldata = new Axis*[col];
        ceils = new CeilEnum*[row];
        
        for (int i = 0; i < row; i++)
            ceils[i] = new CeilEnum[col];
            
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                ceils[i][j] = UNKNOWN;
    }
    
    void setRow(const int r, vector<int> *data) {
        rowdata[r] = new Axis(r, -1, data->size(), col, data);
    }
    
    void setCol(const int c, vector<int> *data) {
        coldata[c] = new Axis(-1, c, data->size(), row, data);
    }
    
    void print(std::ostream &os) const {
        for (int i = 0; i < row; i++) {
            if (trace)
                os << std::setw(2) << i << ": ";
            rowdata[i]->display(os);
            if (trace)
                rowdata[i]->print(os);
            os << endl;
        }
    }
    
    bool checkIfSolved() {
        for (int c = 0; c < col; c++) {
            CeilEnum *ceilsData = coldata[c]->getCeils();

            for (int r = 0; r < row; r++)
                ceilsData[r] = ceils[r][c];

            if (! coldata[c]->checkByCeils())
                return false;
        }
        return true;
    }

    bool dfs_traverse (const int r) {

#if ENABLE_CTRL_C_TRACE
        if (trace) {
            cerr << "Trace on row: " << r << endl;
            this->print(cerr);
            cerr << endl;
            trace = false;
        }
#endif
        
        if (r == row) {
            return checkIfSolved();
        } else {
            // rowdata[r]->pushState();
            
            CeilEnum *combination = allowed_combination (r);
            
            rowdata[r]->pushState();
        
            while (rowdata[r]->nextState()) {
                // if (trace) {
                //     cerr << "State on row " << r << " " << *rowdata[r]->getState() << endl;
                // }
                if (rewrite_row(r) && early_checker (r, combination)) {
                    if (dfs_traverse(r + 1))
                        return true;
                }
            }
        
            delete combination;
        
            rowdata[r]->popState();
            return false;
        }
    }
    
    bool solve() {
        initial_guess();

        return dfs_traverse (0);
    }
};

std::ostream& operator << (std::ostream& os, Nonogram& nonogram) {
    nonogram.print(os);
    return os;
}

Nonogram* Nonogram::inputReader () {
    Nonogram *nonogram = NULL;
    string line;
    
    if (std::getline (std::cin, line)) {
        vector<int> *list = Utilities::toIntArray(line);
        int row = list->at(0);
        int col = list->at(1);
        list->clear();
        delete list;
        
        nonogram = new Nonogram(row, col);
        
        for (int r = 0; r < row; r++) {
            std::getline(std::cin, line);
            nonogram->setRow(r, Utilities::toIntArray(line));
        }
        
        for (int c = 0; c < col; c++) {
            std::getline(std::cin, line);
            nonogram->setCol(c, Utilities::toIntArray(line));
        }
    }

    return nonogram;
}

int main (int argc, char *argv[]) {
    Nonogram *nonogram = Nonogram::inputReader ();
    
#if ENABLE_CTRL_C_TRACE
    signal (SIGINT, ctrlc_handler);
#endif

    bool solved = nonogram->solve();
    
    if (solved)
        cout << *nonogram;

    delete nonogram;

    return solved ? 0 : -1;
}
