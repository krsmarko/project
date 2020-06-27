/*! \file       problem_a.cpp
 *  \author     Student Name
 *  \version    0.1
 *  \date       15.01.2019
 *
 *  Solution of the Problem A/Contest 1.
 *
 *  Task: <paste it here if you need it>
 *
 *  Input format: <paste it here if you need it>
 *
 *  Output format: <paste it here if you need it>
 */


#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ctime>
#include <unordered_set>
#include <queue>
#include <unistd.h>

int counter = 0;


class widthEl
{
public:
    widthEl(int num, int prev):
        num_(num),
        prev_(prev)
    {
    }



    int num_;
    int prev_;

    ~widthEl()
    {}


};

class Term
{

public:
    Term(int  number, int deg, double coef):
            number_(number),
            deg_(deg),
            coef_(coef)
    {

    }


    int number_;
    int deg_;
    double coef_;
    std::unordered_map<int, Term*> umap_;



};

class HummerFunc
{
public:
    bool checkTerm (std::vector<int>& term);
    Term* addTerm (std::vector<int> term, double coef = 0);
    Term* getTerm (std::vector<int>& term);
    double calculate (std::vector<int>& input, Term* start);
    void khumavala1 (std::vector<int>& input, Term* start, std::vector<double>& output, std::vector<int>& deep);
    double lowBound (std::vector<int>& input, Term* start, std::vector<int>& deep);
    double lowBoundStrong (std::vector<int>& input, Term* start, std::vector<int>& deep,
                           std::vector<double>& khumavala2);
    double  elenkotterLowBound (std::vector<int>& input);
    void khumavala2 (std::vector<int>& input, Term* start, std::vector<double>& output, std::vector<int>& deep);
    void dataCorrecting (std::vector<int>& input, Term* start, std::vector<int>& deep, double upperBound, int* prev);
    void zeroTermInfo (std::vector<int>& input, Term* start, std::vector<int>& deep, std::vector<int>& zero, std::vector<int>& nonZero);
    void MQ (std::vector<int>& input, Term* start, std::vector<std::vector<double>>& outputQCoef, std::vector<std::vector<double>>& outputQIN, std::vector<int>& deep);
    void phantomKhumavala (std::vector<int>& input, Term* start, std::vector<double>& khumavala1,
                           std::vector<double>& khumavala2, std::vector<int>& deep,
                           std::vector<int>& phantomDeep0, std::vector<int>& phantomDeep1);
    void phantomKhumavala1 (std::vector<int>& input, Term* start, std::vector<double>& khumavala1,
                            std::vector<int>& deep);

    //double calcExtT (Term& start, int key, std::vector<int>& input);
    void printHF (std::vector<int>& input, Term* start, std::vector<int>& deep);
    void opt(std::vector<int> input, std::vector<int>& processData, double& minValue, std::vector<int>& answer);
    bool inputIsFull (std::vector<int>& input);
    bool inputIsUnity (std::vector<int>& input);

    HummerFunc ():
            zero(-1, 0, 0)
    {}


    Term zero;
    std::map<std::vector<int>, Term*> terms;
    std::unordered_set<int> was;
    std::unordered_set<int> wasC;

    ~HummerFunc()
    {
        for(auto i: terms)
        {
           // delete i.second;
        }
    }
};

int main()
{
    HummerFunc func;


    std::ifstream fin("cap131.txt");
    int m, n;
    fin >> m >> n;
    int trash;


    std::vector<double> costs;

    for(int i = 0; i < m; i++)
    {
        double tempCoef;
        fin >> trash;
        fin >> tempCoef;
        costs.push_back(tempCoef);
    }

    std::vector<std::vector<std::pair<double, int>>> data;

    for(int i = 0; i < n; i++)
    {
        std::vector<std::pair<double, int>> tempColumn;
        data.push_back(tempColumn);
        fin >> trash;
        for(int j = 0; j < m; j++)
        {
            double tempCoef;
            fin >> tempCoef;
            data[i].push_back(std::pair<double, int> (tempCoef, j));
        }
    }



    /*

    std::vector<double> costs;

    std::vector<std::vector<std::pair<double, int>>> data;

    for(int i = 0; i < n; i++)
    {
        std::vector<std::pair<double, int>> tempColumn;
        data.push_back(tempColumn);
    }

    for(int i = 0; i < m; i++)
    {
        double tempCoef;
        fin >> trash;
        fin >> tempCoef;
        costs.push_back(tempCoef);

        std::vector<std::pair<double, int>> tempColumn;
        data.push_back(tempColumn);
        for(int j = 0; j < n; j++)
        {
            double tempCoef;
            fin >> tempCoef;
            data[j].push_back(std::pair<double, int> (tempCoef, i));
        }
    }
     */



    fin.close();

    clock_t start = std::clock();

    std::vector<int> temp = {0};

    for(int i = 0; i < m; i++)
    {
        temp[0] = i;
        func.zero.coef_ += costs[i];
        func.addTerm(temp, -costs[i]);
    }

    for(int i = 0; i < n; i++)
    {
        std::sort(data[i].begin(), data[i].end());


        func.zero.coef_ += data[i][0].first;
        temp = {};
        temp.push_back(data[i][0].second);
        for(int j = 1; j < m; j++)
        {

            func.addTerm(temp, data[i][j].first - data[i][j - 1].first);
            temp.push_back(data[i][j].second);


        }

    }

    int min = 0;
    int max = 0;
    int sum = 0;

    int counter = 0;
    /*
    for(auto i = func.terms.begin(); i != func.terms.end(); i++)
    {

        if(min > i->second->coef_)
        {
            min = i->second->coef_;
        }

        if(max < i->second->coef_)
        {
            max = i->second->coef_;
        }

        if(100 <= i->first.size())
        {
            counter++;
            sum += i->second->coef_;
        }


    }

    std::cout << max << " " << min << " " << counter << " " << double(sum) << "\n";
     */

    costs.clear();

    data.clear();






    /*
    std::vector<double> output = {0, 0, 0, 0};
    std::vector<int> deep;

    std::vector<int> term = {0, 3};
    //func.printHF(temp, 66, &func.zero, deep);

    func.khumavala2(temp, 2, &func.zero, output, deep);


    for(int i = 0; i < 4; i++)
    {
        std::cout << " " << output[i] << " ";
    }

    temp = {1, 0, 0, 1};

    */
    //std::cout << " " << func.calculate(temp, 1, &func.zero);


    std::vector<int> input(m, -1);
    std::vector<int> input0(m, 0);
    std::vector<int> answerVec(m, 0);
    std::vector<int> processData(4, 0);
    std::vector<double> output(m, 0);
    std::vector<int> deep;

    std::vector<double> output1(input.size(), 0.0);
    std::vector<double> output2(input.size(), 0.0);

    int countK1 = 0;
    int countK2 = 0;





    func.khumavala1(input, &func.zero, output1, deep);

    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == -1 && output1[i] >= 0)
        {
            countK1++;
            input[i] = 0;
        }
    }


    func.khumavala2(input, &func.zero, output2, deep);


    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == -1 && output2[i] <= 0)
        {
            countK2++;
            input[i] = 1;
        }
    }





    int extraK2 = countK2;
    int extraK1 = countK1;

    while(extraK2 > 0)
    {

        std::fill(output1.begin(), output1.end(), 0.0);
        std::fill(output2.begin(), output2.end(), 0.0);
        extraK2 = 0;
        extraK1 = 0;
        func.khumavala1(input, &func.zero, output1, deep);



        for(int i = 0; i < input.size(); i++)
        {
            if(input[i] == -1 && output1[i] >= 0)
            {
                extraK1++;
                countK1++;
                input[i] = 0;
            }
        }

        if(extraK1 == 0)
            break;

        func.khumavala2(input, &func.zero, output2, deep);


        for(int i = 0; i < input.size(); i++)
        {
            if(input[i] == -1 && output2[i] <= 0)
            {
                countK2++;
                extraK2++;
                input[i] = 1;
            }
        }
    }

    func.inputIsFull(input);

    processData[2] += countK1;
    processData[3] += countK2;

    //std::cout << "\n b: " << deep.size();




   // std::cout << "init1 " << func.elenkotterLowBound(input) << "\n";


    std::fill(output1.begin(), output1.end(), 0.0);
    std::fill(output2.begin(), output2.end(), 0.0);


    func.khumavala1(input, &func.zero, output1, deep);
    func.khumavala2(input, &func.zero, output2, deep);

    //std::cout << "\n" << func.elenkotterLowBound(input) << '\n';

    std::vector<int> zero1(input.size(), 0.0);
    std::vector<int> nonZero1(input.size(), 0.0);

    //func.zeroTermInfo (input, &func.zero, deep, zero1, nonZero1);
    /*
    for(int i = 0; i < m; i++)
    {
        std::cout << i + 1 << " " << zero1[i] << " " << nonZero1[i] << "\n";
    }
     */

    int prev = -1;
    func.dataCorrecting(input, &func.zero, deep, 793439.5625, &prev);

    //std::cout << "\n" << func.elenkotterLowBound(input) << '\n';

    //std::cout << "\n c: " << deep.size();



    std::vector<double> output11(input.size(), 0.0);
    std::vector<double> output22(input.size(), 0.0);

    func.khumavala1(input, &func.zero, output11, deep);
    func.khumavala2(input, &func.zero, output22, deep);


    std::vector<int> zero(input.size(), 0.0);
    std::vector<int> nonZero(input.size(), 0.0);

    func.zeroTermInfo (input, &func.zero, deep, zero, nonZero);
    /*
    for(int i = 0; i < m; i++)
    {
        std::cout << i + 1 << " " << zero[i] << " " << nonZero[i] << "\n";
    }
     */





    /*

    for(int i = 0; i < m; i++)
    {
        std::cout << output1[i] << " " << output11[i] << " " << output2[i] << " " << output22[i] <<"\n";
    }
     */

    //std::cout.precision(13);

    //std::cout << "\n" << func.zero.coef_ << "\n";

    //func.printHF(input, &func.zero, deep);

    for(int i = 0; i < m; i++)
    {
        if(input[i] == -1)
        {
            answerVec[i] = 0;
        }
        else
        {
            answerVec[i] = input[i];
        }
    }



    double answer = func.calculate(answerVec, &func.zero);
    //std::cout << "init" << answer << " " << func.elenkotterLowBound(input) << "\n";
    func.inputIsFull(input);
    /*
    for(int i = 0; i < 30; i++)
    {
        input[i] = 1;
    }
     */



    func.opt(input, processData, answer, answerVec);
    /*
    std::fill(output1.begin(), output1.end(), 0.0);
    std::fill(output2.begin(), output2.end(), 0.0);


    func.khumavala1(input, &func.zero, output1, deep);
    func.khumavala2(input, &func.zero, output2, deep);







    for(int i = 0; i < m; i++)
    {
        std::cout << output1[i] << " " << output2[i] << "\n";
    }
    */




    if(func.inputIsUnity(answerVec))
    {
        func.khumavala2(input, &func.zero, output, deep);

        int futureZero = -1;

        for(int i = 0; i < input.size(); i++)
        {
            if(futureZero == -1 || output[futureZero] < output[i])
                futureZero = i;
        }

        answerVec[futureZero] = 0;
    }





    clock_t end = std::clock();

    std::cout.precision(13);
    func.inputIsFull(answerVec);

    std::cout << double(answer) << "\n" << func.calculate(answerVec, &func.zero) <<  "\n";


    for(int i = 0; i < 4; i++)
    {
        std::cout << processData[i] << "\n";
    }

    std::cout << end - start << "\n";
    std::cout << (end - start)/float(CLOCKS_PER_SEC) << "\n";

    for(int i = 0; i < m; i++)
    {
        std::cout << answerVec[i] << " ";
    }





}

bool HummerFunc::checkTerm (std::vector<int>& term)
{
    Term& currentTerm = zero;

    for(int i = 0; i < term.size(); i++)
    {
        if(currentTerm.umap_.find(term[i]) == currentTerm.umap_.end())
        {
            return false;
        }

        currentTerm = *(currentTerm.umap_.find(term[i])->second);

    }

    return true;
}

Term* HummerFunc::addTerm (std::vector<int> term, double coef)
{
    Term* currentTerm = &zero;


    for(int i = 0; i < term.size() - 1; i++)
    {

        currentTerm = currentTerm->umap_.find(term[i])->second;
    }


    int temp = *(term.end() - 1);

    std::sort(term.begin(), term.end());

    if(terms.find(term) == terms.end())
    {
        Term* newTerm = new Term(terms.size(), term.size(), 0);
        terms[term] = newTerm;

    }


    currentTerm->umap_[temp] = terms[term];
    currentTerm->umap_[temp]->coef_ += coef;

    return terms[term];
}

Term* HummerFunc::getTerm (std::vector<int>& term)
{
    Term* currentTerm = &zero;

    for(int i = 0; i < term.size(); i++)
    {
        currentTerm = currentTerm->umap_.find(term[i])->second;
    }

    return currentTerm;
}

double HummerFunc::calculate (std::vector<int>& input, Term* start)
{
    if(start->number_ == -1)
        was.clear();

    double answer = 0;

    if(was.find(start->number_) == was.end())
    {
        for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
        {
            if(input[i->first] != 0)
                answer += calculate(input, i->second);
        }

        answer += start->coef_;

        was.insert(start->number_);
    }


    return answer;
}



void HummerFunc::khumavala1 (std::vector<int>& input, Term* start, std::vector<double>& output, std::vector<int>& deep)
{
    if(start->number_ == -1)
        was.clear();

    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(was.find(i->second->number_) == was.end())
        {

            if (input[i->first] == -1)
            {
                was.insert(i->second->number_);

                if(!deep.empty())
                    return;

                deep.push_back(i->first);


                output[deep[0]] += i->second->coef_;

                khumavala1(input, i->second, output, deep);


                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                khumavala1(input, i->second, output, deep);

                if(deep.size() == 1)
                {
                    output[deep[0]] += i->second->coef_;


                }

                was.insert(i->second->number_);
            }
        }
    }

}

void HummerFunc::khumavala2 (std::vector<int>& input, Term* start, std::vector<double>& output, std::vector<int>& deep)
{


    if(start->number_ == -1)
        was.clear();

    //std::cout << "size: " << deep.size() << "\n";

    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(was.find(i->second->number_) == was.end())
        {

            if (input[i->first] == -1)
            {
                was.insert(i->second->number_);

                deep.push_back(i->first);

                for(int j = 0; j < deep.size(); j++)
                {
                    output[deep[j]] += i->second->coef_;

                }

                //std::cout << i->second->coef_ << "\n";

                khumavala2(input, i->second, output, deep);


                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                khumavala2(input, i->second, output, deep);

                if(!deep.empty())
                {
                    for(int j = 0; j < deep.size(); j++)
                    {
                        output[deep[j]] += i->second->coef_;

                    }

                    //std::cout << i->second->coef_ << "\n";


                }

                was.insert(i->second->number_);
            }
        }
    }


}

void HummerFunc::phantomKhumavala1 (std::vector<int>& input, Term* start, std::vector<double>& khumavala1,
                                    std::vector<int>& deep)
{

    if(start->number_ == -1)
        was.clear();

    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(was.find(i->second->number_) == was.end())
        {


            if (input[i->first] == 1 )
            {
                deep.push_back(i->first);

                phantomKhumavala1 (input, i->second, khumavala1, deep);



                for(int j = 0; j < deep.size(); j++)
                {
                    khumavala1[deep[j]] += i->second->coef_;
                }


                was.insert(i->second->number_);

                deep.pop_back();

            }

        }
    }
}


void HummerFunc::phantomKhumavala (std::vector<int>& input, Term* start, std::vector<double>& khumavala1,
                                   std::vector<double>& khumavala2, std::vector<int>& deep,
                                   std::vector<int>& phantomDeep0, std::vector<int>& phantomDeep1)
{

    if(start->number_ == -1)
        was.clear();

    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(was.find(i->second->number_) == was.end())
        {

            if (input[i->first] == -1)
            {
                was.insert(i->second->number_);

                deep.push_back(i->first);

                if(!phantomDeep0.empty())
                {
                    khumavala2[phantomDeep0[0]] += i->second->coef_;
                }

                phantomKhumavala (input, i->second, khumavala1,
                                  khumavala2, deep,
                                  phantomDeep0, phantomDeep1);


                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                phantomDeep1.push_back(i->first);

                phantomKhumavala (input, i->second, khumavala1,
                                  khumavala2, deep,
                                  phantomDeep0, phantomDeep1);


                if(!phantomDeep0.empty())
                {
                    khumavala2[phantomDeep0[0]] += i->second->coef_;
                }
                else if(deep.empty())
                {
                    for(int j = 0; j < phantomDeep1.size(); j++)
                    {
                        khumavala1[phantomDeep1[j]] += i->second->coef_;
                    }

                }

                was.insert(i->second->number_);

                phantomDeep1.pop_back();

            }
            else if (input[i->first] == 0 && phantomDeep0.empty())
            {

                phantomDeep0.push_back(i->first);

                phantomKhumavala (input, i->second, khumavala1,
                                  khumavala2, deep,
                                  phantomDeep0, phantomDeep1);



                khumavala2[phantomDeep0[0]] += i->second->coef_;


                was.insert(i->second->number_);

                phantomDeep0.pop_back();

            }

        }
    }
}



//double HummerFunc::calcExtT (Term& start, int key, std::vector<int>& input);
void HummerFunc::printHF (std::vector<int>& input, Term* start, std::vector<int>& deep)
{
    if(start->number_ == -1)
    {
        was.clear();
        std::cout << start->coef_;
    }

    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(was.find(i->second->number_) == was.end())
        {

            if (input[i->first] == -1)
            {
                deep.push_back(i->first);

                std::cout << "+" << i->second->coef_;
                for(int j = 0; j < deep.size(); j++)
                {
                    std::cout<< "y" << deep[j] + 1;
                }
                printHF(input, i->second, deep);

                was.insert(i->second->number_);


                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                printHF(input, i->second, deep);


                    std::cout << "+" << i->second->coef_;
                    for(int j = 0; j < deep.size(); j++)
                    {
                        std::cout<< "y" << deep[j] + 1;
                    }




                was.insert(i->second->number_);
            }
        }
    }


}

double HummerFunc::lowBound (std::vector<int>& input, Term* start, std::vector<int>& deep)
{
    double answer = 0;

    if(start->number_ == -1)
    {
        answer += start->coef_;

        was.clear();
    }




    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(was.find(i->second->number_) == was.end())
        {

            if (input[i->first] == -1)
            {
                was.insert(i->second->number_);

                if(!deep.empty())
                    return answer;

                deep.push_back(i->first);


                answer += i->second->coef_;

                answer += lowBound(input, i->second, deep);


                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                answer += lowBound(input, i->second, deep);

                if(deep.size() <= 1)
                {
                    answer += i->second->coef_;


                }

                was.insert(i->second->number_);
            }
        }
    }

    return answer;

}

double HummerFunc::lowBoundStrong (std::vector<int>& input, Term* start, std::vector<int>& deep,
                                   std::vector<double>& khumavala2)
{
    double answer = 0;

    if(start->number_ == -1)
    {
        answer += start->coef_;

        was.clear();
    }




    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(was.find(i->second->number_) == was.end())
        {

            if (input[i->first] == -1)
            {
                was.insert(i->second->number_);

                deep.push_back(i->first);


                answer += i->second->coef_;

                answer += lowBoundStrong(input, i->second, deep, khumavala2);


                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                answer += lowBoundStrong(input, i->second, deep, khumavala2);


                answer += i->second->coef_;




                was.insert(i->second->number_);
            }
        }



    }


    if(deep.size() > 1)
    {
        double maxKhumavala = 0;

        for (int j = 0; j < deep.size(); j++) {
            if (khumavala2[deep[j]] > maxKhumavala) {
                maxKhumavala = khumavala2[deep[j]];
            }

            khumavala2[deep[j]] -= start->coef_;
        }

        answer -= std::min(maxKhumavala, start->coef_);
    }

    return answer;

}

void HummerFunc::zeroTermInfo (std::vector<int>& input, Term* start, std::vector<int>& deep, std::vector<int>& zero, std::vector<int>& nonZero)
{
    if(start->number_ == -1)
        wasC.clear();


    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(wasC.find(i->second->number_) == wasC.end())
        {

            if (input[i->first] == -1)
            {

                wasC.insert(i->second->number_);


                deep.push_back(i->first);



                zeroTermInfo(input, i->second, deep, zero, nonZero);

                if(deep.size() > 1)
                {
                    if(i->second->coef_ == 0)
                        zero[deep.size()]++;
                    else
                        nonZero[deep.size()]++;
                }






                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                wasC.insert(i->second->number_);
                zeroTermInfo(input, i->second, deep, zero, nonZero);

                if(deep.size() > 1)
                {
                    if(i->second->coef_ == 0)
                        zero[deep.size()]++;
                    else
                        nonZero[deep.size()]++;
                }



            }
        }
    }


}


void HummerFunc::dataCorrecting (std::vector<int>& input, Term* start, std::vector<int>& deep, double upperBound, int* prev)
{
    if(start->number_ == -1)
        wasC.clear();

    //std::cout << "Start:" << deep.size() << "\n";

    //inputIsFull(input);

    //std::cout << "size: " << deep.size() << "\n";

    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(wasC.find(i->second->number_) == wasC.end())
        {

            if (input[i->first] == -1)
            {

                wasC.insert(i->second->number_);


                deep.push_back(i->first);



                dataCorrecting(input, i->second, deep, upperBound, prev);






                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                wasC.insert(i->second->number_);
                dataCorrecting(input, i->second, deep, upperBound, prev);



            }
        }
    }



    //std::cout << "end:" << deep.size() << "\n";

    if(deep.size() > 1)
    {
        for(int i = 0; i < deep.size(); i++)
        {
            input[deep[i]] = 1;
        }

        if(*prev != -1)
            //input[*prev] = 0;




        int countK1 = 0;
        int countK2 = 0;

        std::vector<double> output1(input.size(), 0.0);
        std::vector<double> output2(input.size(), 0.0);
        std::vector<int> deepK;
        std::vector<int> correctedVariable;
        /*
        khumavala1(input, &zero, output1, deepK);

        for (int i = 0; i < input.size(); i++)
        {
            if (input[i] == -1 && output1[i] >= 0) {
                countK1++;
                input[i] = 0;
                correctedVariable.push_back(i);
            }
        }


        khumavala2(input, &zero, output2, deepK);


        for (int i = 0; i < input.size(); i++)
        {
            if (input[i] == -1 && output2[i] <= 0) {
                countK2++;
                input[i] = 1;
                correctedVariable.push_back(i);
            }
        }


        int extraK2 = countK2;
        int extraK1 = countK1;

        while (extraK2 > 0) {
            std::fill(output1.begin(), output1.end(), 0.0);
            std::fill(output2.begin(), output2.end(), 0.0);
            extraK2 = 0;
            extraK1 = 0;
            khumavala1(input, &zero, output1, deepK);


            for (int i = 0; i < input.size(); i++) {
                if (input[i] == -1 && output1[i] >= 0) {
                    extraK1++;
                    countK1++;
                    input[i] = 0;
                    correctedVariable.push_back(i);
                }
            }

            if (extraK1 == 0)
                break;

            khumavala2(input, &zero, output2, deepK);


            for (int i = 0; i < input.size(); i++) {
                if (input[i] == -1 && output2[i] <= 0) {
                    countK2++;
                    extraK2++;
                    input[i] = 1;
                    correctedVariable.push_back(i);
                }
            }
        }



        if (extraK2 > 0)
        {
            std::fill(output2.begin(), output2.end(), 0.0);
            khumavala2(input, &zero, output2, deepK);
        }
         */






        //double lowBound = lowBoundStrong(input, &zero, deepK, output2);

       // std::cout << "lowBound:" << lowBound << "\n";
        //double lowBoundl = lowBound(input, &zero, deepK);
        double lowBound = elenkotterLowBound(input);
        std::fill(output1.begin(), output1.end(), 0.0);


        phantomKhumavala1(input, &zero, output1, deepK);
        /*
        double max = output1[deep[0]];

        for(int i = 0; i < deep.size(); i++)
        {
            if(max < output1[deep[i]])
            {
                max = output1[deep[i]];
            }

        }
         */


        double corectCoef = lowBound - upperBound;

        //printf("%f\n", corectCoef - max);


        corectCoef = corectCoef - 0.01;

        if(start->coef_ < 0)
        {
            std::cout << "rrrrrrrrrrrrrrrrrrrrrrrrrrrr";
        }


        if (corectCoef > 0)
        {

            if(std::min(corectCoef, start->coef_) < 0)
            {
                std::cout << "rrrrrrrrrrrrrrrrrrrrrrrrrrrr";
            }

            //std::cout << "size: " << deep.size();
            double change = std::min(corectCoef, start->coef_);
            start->coef_ -= std::min(corectCoef, start->coef_);

            corectCoef = elenkotterLowBound(input) - upperBound;


                //std::cout << "\n" << elenkotterLowBound(input) << " " << change << " " << lowBound << "\n";
                //printf("\n%f %f %f \n", elenkotterLowBound(input), change, lowBound);



        }


        //printf("\n%f %f %f\n", elenkotterLowBound(input), lowBound, upperBound);




        for(int i = 0; i < deep.size(); i++)
        {
            input[deep[i]] = -1;
        }

        //printf("\n%f %f \n", elenkotterLowBound(input), lowBound);

        for (int i = 0; i < correctedVariable.size(); i++)
        {
            input[correctedVariable[i]] = -1;
        }

        if(*prev != -1)
            input[*prev] = -1;

        if(start->coef_ == 0 && corectCoef != 0)
        {
            *prev = deep.back();
        }


        //printf("\n%f %f \n", elenkotterLowBound(input), lowBound);
        counter++;
        std::cout << counter << "\n";




    }





}


void HummerFunc::MQ (std::vector<int>& input, Term* start, std::vector<std::vector<double>>& outputQCoef, std::vector<std::vector<double>>& outputQIN, std::vector<int>& deep)
{
    if(start->number_ == -1)
        was.clear();

    for(auto i = start->umap_.begin(); i != start->umap_.end(); i++)
    {


        if(was.find(i->second->number_) == was.end())
        {

            if (input[i->first] == -1)
            {
                was.insert(i->second->number_);


                deep.push_back(i->first);

                if(deep.size() == 2)
                {

                    outputQCoef[deep[0]][deep[1]] += i->second->coef_;
                    outputQCoef[deep[1]][deep[0]] += i->second->coef_;

                    outputQIN[deep[0]][deep[1]] += i->second->coef_;
                    outputQIN[deep[1]][deep[0]] += i->second->coef_;
                }
                else if(deep.size() > 2)
                {
                    for (int c = 0; c < deep.size(); c++)
                    {
                        for (int j = 0; j < deep.size(); j++)
                        {
                            outputQIN[deep[c]][deep[j]] += i->second->coef_;
                        }
                    }

                }

                MQ(input, i->second, outputQCoef, outputQIN, deep);


                deep.pop_back();


            }
            else if (input[i->first] == 1 )
            {
                MQ(input, i->second, outputQCoef, outputQIN, deep);

                if(deep.size() == 2)
                {

                    outputQCoef[deep[0]][deep[1]] += i->second->coef_;
                    outputQCoef[deep[1]][deep[0]] += i->second->coef_;

                    outputQIN[deep[0]][deep[1]] += i->second->coef_;
                    outputQIN[deep[1]][deep[0]] += i->second->coef_;


                }
                else if(deep.size() > 2)
                {
                    for (int c = 0; c < deep.size(); c++)
                    {
                        for (int j = 0; j < deep.size(); j++)
                        {
                            outputQIN[deep[c]][deep[j]] += i->second->coef_;
                        }
                    }


                }
                was.insert(i->second->number_);

            }
        }
    }

}

void HummerFunc::opt(std::vector<int> input, std::vector<int>& processData, double& minValue, std::vector<int>& answer)
{


    if(inputIsFull(input))
    {

        processData[1]++;

        double tempValue = calculate(input, &zero);


        if(tempValue < minValue)
        {

            minValue = tempValue;
            std::cout << tempValue << "\n";
            answer = input;
        }

        return;

    }


    int countK1 = 0;
    int countK2 = 0;

    std::vector<double> output1(input.size(), 0.0);
    std::vector<double> output2(input.size(), 0.0);
    std::vector<int> deep;

    khumavala1(input, &zero, output1, deep);

    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == -1 && output1[i] >= 0)
        {
            countK1++;
            input[i] = 0;
        }
    }


    khumavala2(input, &zero, output2, deep);


    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == -1 && output2[i] <= 0)
        {
            countK2++;
            input[i] = 1;
        }
    }





    int extraK2 = countK2;
    int extraK1 = countK1;

    while(extraK2 > 0)
    {

        std::fill(output1.begin(), output1.end(), 0.0);
        std::fill(output2.begin(), output2.end(), 0.0);
        extraK2 = 0;
        extraK1 = 0;
        khumavala1(input, &zero, output1, deep);



        for(int i = 0; i < input.size(); i++)
        {
            if(input[i] == -1 && output1[i] >= 0)
            {
                extraK1++;
                countK1++;
                input[i] = 0;
            }
        }

        if(extraK1 == 0)
            break;

        khumavala2(input, &zero, output2, deep);


        for(int i = 0; i < input.size(); i++)
        {
            if(input[i] == -1 && output2[i] <= 0)
            {
                countK2++;
                extraK2++;
                input[i] = 1;
            }
        }
    }



    processData[2] += countK1;
    processData[3] += countK2;

    if(countK1 > 0 || countK2 > 0)
    {
        if (inputIsFull(input))
        {
            processData[1]++;


            double tempValue = calculate(input, &zero);


            if(tempValue < minValue)
            {
                minValue = tempValue;
                std::cout << tempValue << "\n";
                answer = input;
            }
            return;
        }

    }



    if(extraK2 > 0)
    {
        std::fill(output2.begin(), output2.end(), 0.0);
        khumavala2(input, &zero, output2, deep);
    }

    //double temp = lowBoundStrong(input, &zero, deep, output2);
    //double temp = lowBound(input, &zero, deep);

    /*
    std::cout << "\n";
    printHF(input, &zero, deep);
    std::cout << "\n";

    */

    double temp = elenkotterLowBound(input);



    //std::cout << "B:" << temp << " " << lowBoundStrong(input, &zero, deep, output2) << "\n";
    if(temp >= minValue)
    {
       // std::cout << "bound";
        //std::cout << "cut";
        return;
    }



    if(extraK1 > 0)
    {
        std::fill(output1.begin(), output1.end(), 0.0);
        khumavala1(input, &zero, output1, deep);
    }





    int next = -1;



    /*
    std::vector<std::vector<double>> output3;


    for(int i = 0; i < input.size(); i++)
    {
        output3.push_back(std::vector<double> (input.size(), 0.0));
    }
     */

    std::vector<double> output4(input.size(), 0.0);
    std::vector<double> output5(input.size(), 0.0);
    std::vector<int> phantomDeep0;
    std::vector<int> phantomDeep1;

    phantomKhumavala (input, &zero, output4, output5, deep,
                      phantomDeep0, phantomDeep1);



    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == 0 && output5[i] < 0)
        {
            //std::cout << "derive2 ";
           return;
        }

        if(input[i] == 1 && output4[i] > 0)
        {
            //std::cout << "derive1 ";
           return;
        }


    }


   // std::cout << "GGG";

    /*
    if(extraK1 > 0)
        khumavala1(input, &zero, output1, deep);

    if(extraK2 > 0)
        khumavala2(input, &zero, output2, deep);



    std::vector<std::vector<double>> output3;
    std::vector<std::vector<double>> output6;

    for(int i = 0; i < input.size(); i++)
    {
        output3.push_back(std::vector<double> (input.size(), 0.0));
        output6.push_back(std::vector<double> (input.size(), 0.0));
    }



    //MQ(input, &zero, output3, output6, deep);

    int branchType = -1;
    int nonBinNext = -1;


    int comp;

    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == -1)
        {
            for(int j = 0; j < input.size(); j++)
            {
                if(i != j && input[j] == -1)
                {

                    if (output3[i][j] + output1[j] + output1[i] >= 0)
                    {
                        if(next == -1 || comp < std::max(output2[i], output2[j]))
                        {
                            next = i;
                            nonBinNext = j;
                            comp = std::max(output2[next], output2[nonBinNext]);
                            branchType = 0;
                        }

                    }

                    if (-output6[i][j] + output2[j] + output2[i] <= 0)
                    {
                        if(next == -1 || comp < std::max(-output1[i], -output1[j]))
                        {
                            next = i;
                            nonBinNext = j;
                            comp = std::max(-output1[next], -output1[nonBinNext]);
                            branchType = 1;
                        }

                    }



                }

            }
        }

    }
     */


    /*
    std::vector<std::vector<double>> output3;

    for(int i = 0; i < input.size(); i++)
    {
        output3.push_back(std::vector<double> (input.size(), 0.0));
    }
    MQ(input, &zero, output3, deep);

    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == -1)
        {
            for(int j = 0; j < input.size(); j++)
            {


                if(i != j && input[j] == -1 && output3[i][j] + output1[j] >= 0)
                {
                    next = i;
                    break;
                }
            }
        }

        if(next != -1)
        {
            //std::cout << "uwu";
            break;
        }
    }
     */

    if(next == -1)
    {
        std::fill(output2.begin(), output2.end(), 0.0);
        khumavala2(input, &zero, output2, deep);


        int comp;


        for (int i = 0; i < input.size(); i++)
        {
            if(next == -1 && input[i] == -1)
            {
                next = i;
                //comp = std::max(output2[next], -output1[next]);
                comp = output2[next];
                //comp = output1[next];
            }

            if (input[i] == -1 && comp < output2[i])
            {
                next = i;
                //comp = std::max(output2[next], -output1[next]);
                comp = output2[next];
                //comp = output1[next];
            }
        }
    }





    output1.clear();
    output2.clear();
    //output3.clear();
    output4.clear();
    output5.clear();
    processData[0] += 2;


    input[next] = 0;
    opt(input, processData, minValue, answer);
    input[next] = 1;
    opt(input, processData, minValue, answer);


}


double  HummerFunc::elenkotterLowBound (std::vector<int>& input)
{

    std::queue<std::pair<Term*, int>> queue;
    std::queue<std::pair<Term*, int>> queuePrior;

    std::vector<widthEl> chain;

    chain.emplace_back(-1, -1);
    std::vector<double> violation (input.size(), 0.0);



    was.clear();

    double answer = zero.coef_;



    queue.emplace(&zero, 0);

    std::pair<Term*, int> current (&zero, 0);


    while(!queue.empty() || !queuePrior.empty())
    {
        //std::cout << "queueu size:"<< queue.size() << " PRqueue:" << queuePrior.size() << "\n";

        if(!queuePrior.empty())
        {
            current.first = queuePrior.front().first;
            current.second = queuePrior.front().second;
            queuePrior.pop();
        }
        else if(!queue.empty())
        {
            current.first = queue.front().first;
            current.second = queue.front().second;
            queue.pop();
        }


        for(auto i = current.first->umap_.begin(); i != current.first->umap_.end(); i++)
        {



            if(was.find(i->second->number_) == was.end())
            {


                if (input[i->first] == -1)
                {

                    was.insert(i->second->number_);



                    chain.emplace_back(i->first, current.second);



                    current.second = chain.size() - 1;



                    if(current.second != 0)
                    {

                        double max = violation[chain[current.second].num_];




                        for(int j = current.second; j != 0;)
                        {


                            if(violation[chain[j].num_] - max > 0)
                            {
                                max = violation[chain[j].num_];
                            }

                            j = chain[j].prev_;
                        }



                        if(max + i->second->coef_ < 0)
                        {
                            for(int j = current.second; j != 0;)
                            {

                                violation[chain[j].num_] += i->second->coef_;
                                j = chain[j].prev_;
                            }

                            answer += i->second->coef_;

                            queue.emplace(i->second, current.second);
                        }
                        else
                        {
                            for(int j = current.second; j != 0;)
                            {
                                violation[chain[j].num_] -= max;
                                j = chain[j].prev_;
                            }

                            answer -= max;
                        }
                    }




                    current.second = chain[current.second].prev_;




                }
                else if (input[i->first] == 1 )
                {


                    was.insert(i->second->number_);


                    if(current.second == 0)
                    {
                        answer += i->second->coef_;

                        queuePrior.emplace(i->second, current.second);
                    }
                    else if(current.second != 0)
                    {
                        double max = violation[chain[current.second].num_];

                        for(int j = current.second; j != 0;)
                        {


                            if(violation[chain[j].num_] - max > 0)
                            {
                                max = violation[chain[j].num_];
                            }

                            j = chain[j].prev_;
                        }

                        if(max + i->second->coef_ < 0)
                        {
                            for(int j = current.second; j != 0;)
                            {

                                violation[chain[j].num_] += i->second->coef_;
                                j = chain[j].prev_;
                            }

                            answer += i->second->coef_;

                            queuePrior.emplace(i->second, current.second);
                        }
                        else
                        {
                            for(int j = current.second; j != 0;)
                            {

                                violation[chain[j].num_] -= max;
                                j = chain[j].prev_;
                            }

                            answer -= max;
                        }
                    }





                }

            }
        }

    }


    return answer;



}


bool HummerFunc::inputIsFull (std::vector<int>& input)
{
    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == -1)
            return false;
    }

    return true;
}


/*
bool HummerFunc::inputIsFull (std::vector<int>& input)
{
    int count = 0;

    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == -1)
            count++;
    }

    std::cout << count << ", ";

    if (count == 0)
        return true;
    else
        return false;
}
*/

bool HummerFunc::inputIsUnity (std::vector<int>& input)
{
    for(int i = 0; i < input.size(); i++)
    {

        if(input[i] != 1)
        {
            return false;
        }
    }

    return true;
}
