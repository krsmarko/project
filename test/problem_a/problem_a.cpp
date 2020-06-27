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
    void khumavala2 (std::vector<int>& input, Term* start, std::vector<double>& output, std::vector<int>& deep);
    void MQ (std::vector<int>& input, Term* start, std::vector<std::vector<double>>& outputQCoef, std::vector<std::vector<double>>& outputQIN, std::vector<int>& deep);
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
    std::vector<int> answerVec(m, 0);
    std::vector<int> processData(4, 0);
    std::vector<double> output(m, 0);
    std::vector<int> deep;

    double answer = func.calculate(answerVec, &func.zero);

    std::cout.precision(13);

    //std::cout << "\n" << func.zero.coef_ << "\n";

    //func.printHF(input, &func.zero, deep);


    func.opt(input, processData, answer, answerVec);




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

    std::cout << double(answer) << "\n" << func.calculate(answerVec, &func.zero) << "\n";


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

                    was.insert(i->second->number_);
                }


            }
        }
    }

}

void HummerFunc::khumavala2 (std::vector<int>& input, Term* start, std::vector<double>& output, std::vector<int>& deep)
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

                for(int j = 0; j < deep.size(); j++)
                {
                    output[deep[j]] += i->second->coef_;
                }

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

                    was.insert(i->second->number_);
                }


            }
        }
    }
}

//double HummerFunc::calcExtT (Term& start, int key, std::vector<int>& input);
void HummerFunc::printHF (std::vector<int>& input, Term* start, std::vector<int>& deep)
{
    if(start->number_ == -1)
        was.clear();

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

                if(!deep.empty())
                {
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

                    was.insert(i->second->number_);
                }


            }
        }
    }

    return answer;

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
        std::fill(output1.begin(), output1.begin(), 0);
        std::fill(output2.begin(), output2.begin(), 0);
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
                answer = input;
            }
            return;
        }

    }


    double temp = lowBound(input, &zero, deep);
    //std::cout << temp << " " << minValue << "\n";
    if(temp >= minValue)
    {
        return;
    }



    if(extraK1 > 0)
        khumavala1(input, &zero, output1, deep);

    if(extraK2 > 0)
        khumavala2(input, &zero, output2, deep);



    std::vector<std::vector<double>> output3;
    std::vector<std::vector<double>> output4;

    for(int i = 0; i < input.size(); i++)
    {
        output3.push_back(std::vector<double> (input.size(), 0.0));
        output4.push_back(std::vector<double> (input.size(), 0.0));
    }



    MQ(input, &zero, output3, output4, deep);

    int branchType = -1;
    int nonBinNext = -1;
    int next = -1;

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
                     
                    if (-output4[i][j] + output2[j] + output2[i] <= 0)
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













    if(next == -1)
    {
        for (int i = 0; i < input.size(); i++) {
            if (next == -1 && input[i] == -1)
                next = i;

            if (input[i] == -1 && std::max(output2[next], -output1[next]) < std::max(output2[i], -output1[i]))
                next = i;
        }
    }




    output1.clear();
    output2.clear();
    output3.clear();
    output4.clear();
    processData[0] += 2;

    if(nonBinNext == -1)
    {
        input[next] = 0;
        opt(input, processData, minValue, answer);
        input[next] = 1;
        opt(input, processData, minValue, answer);
    }
    else
    {
        input[next] = branchType;
        input[nonBinNext] = 1 - branchType;

        opt(input, processData, minValue, answer);

        input[next] = -1;
        input[nonBinNext] = branchType;

        opt(input, processData, minValue, answer);


    }



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