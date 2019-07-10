#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include <math.h>
#include <string.h> 
#include <assert.h>

int interval=8;
int roundnum = 3;
clock_t startclock=clock();
#define po() printf("\nokay\n")
// struct
typedef struct matrix
{
    float** values;
    int row;
    int col;
}matrix;
// genterater function
struct matrix onevalues(int row, int col, float value)
{
    int i=0;
    int j=0;
    matrix result;
    result.row = row;
    result.col = col;
    result.values = (float **)malloc(sizeof(float*) * row);
    for(i=0; i<row; i++)
    {
        result.values[i] = (float *)malloc(sizeof(float) * col);
    }
    for(i=0; i<row; i++)
    {
        for(j=0; j<col; j++)
        {
            result.values[i][j] = value;
        }
    }
    return result;
}
// Utills
void print_uniformally(float cont)
{
    if(cont==0)
    {
        printf(" 0%*s", interval-1, "");
    }
    else if(cont > 0)
    {
        printf(" %-*.*f ", interval-1, interval-5, cont);
    }
    else
    {
        printf("%-*.*f ", interval, interval-5, cont);
    }
}
void print_matrix(matrix foo)
{
    printf("\n");
    for(int i=0; i<foo.row; i++)
    {
        printf("|");
        for(int j=0; j<foo.col; j++)
        {
            print_uniformally(foo.values[i][j]);
        }
        printf("|");
        printf("\n");
    }
    printf("\n");
}
void print_equation(matrix coeffi, matrix consta)
{
    printf("\n");
    for(int i=0; i<coeffi.row; i++)
    {
        printf("|");
        for(int j=0; j<coeffi.col+1; j++)
        {
            if(j==coeffi.col)
            {
                printf("| |");
                print_uniformally(consta.values[i][0]);
            }
            else
            {
                print_uniformally(coeffi.values[i][j]);
            }
        }
        printf("|");
        printf("\n");
    }
}
struct matrix copy_matrix(matrix mat)
{
    matrix result;
    result = onevalues(mat.row, mat.col, 0.0);
    for(int i=0; i<mat.row; i++)
    {
        for(int j=0; j<mat.col; j++)
        {
            result.values[i][j] = mat.values[i][j];
        }
    }
    return result;
}
struct matrix free_matrix(matrix mat)
{
    for(int i=0; i<mat.row; i++)
    {
        free(mat.values[i]);
    }
    free(mat.values);
}
struct matrix input_matrix()
{
    int row = 0;
    int col = 0;
    int ii=0;
    char input_array[2*col];
    matrix result;
    printf("\nInput number of rows : ");
    scanf("%d", &row);
    printf("\nInput number of columns : ");
    scanf("%d", &col);
    result = onevalues(row, col, 0.0);
    for(int i=0; i<row; i++)
    {
        printf("\nInput %dth row of matrix : \n", i+1);
        for(int j=0; j<col; j++)
        {
            scanf("%f", &result.values[i][j]);
        }
    }
    return result;
}
struct matrix clip_by_col(matrix mat, int i, int j)
{
    if(i<0 or i>mat.col or j<0 or j>mat.col or i>=j)
    {
        printf("\nwrong index\n");
        return mat;
    }
    int row = mat.row;
    int col = j-i;
    matrix result = onevalues(row, col, 0.0);
    for(int jj=0; jj<col; jj++)
    {
        for(int ii=0; ii<row; ii++) result.values[ii][jj] = mat.values[ii][jj+i];
    }
    return result;
}
float msur_time(bool startq)
{
    if(startq) startclock = clock();
    else return (float)(clock() - startclock)/CLOCKS_PER_SEC;
    return 0.0;
}
// constant operation
int randint(int min, int max)
{
    int random = rand();
    random = (random % (max-min)) + min;
    return random;
}
float n_mul(float i, int n)
{
    double result=1.0;
    for(int j=0; j<n; j++)
    {
        result *= i;
    }
    return result;
}
float roundpint(float num)
{
    float temp = n_mul(10.0, roundnum);
    return roundf(num*temp)/temp;
}
struct matrix tidy_matrix(matrix mat)
{
    float temp;
    matrix result=copy_matrix(mat);
    for(int i=0; i<mat.row; i++)
    {
        for(int j=0; j<mat.col; j++)
        {
            temp = n_mul(10.0, roundnum);
            result.values[i][j] = roundf(result.values[i][j]*temp)/temp;
        }
    }
    return result;
}
// replacer
struct matrix replace_column(matrix mat1, matrix column, int j)
{
    matrix result;
    result = copy_matrix(mat1);
    for(int i=0; i<mat1.row; i++)
    {
        result.values[i][j] = column.values[i][0];
    }
    return result;
}
//checker
bool check_zerorow(matrix mat)
{
    for(int i=0; i<mat.row; i++)
    {
        for(int j=0; j<mat.col; j++){
            if(mat.values[i][j]!=0)
            {
                break;
            }
            if(j==mat.col-1)
            {
                return true;
            }
        }
    }
    return false;
}
bool check_identity(matrix mat)
{
    if(mat.row!=mat.col) return false;
    for(int i=0; i<mat.row; i++)
        for(int j=0; j<mat.col; j++)
            if(!((i==j and mat.values[i][j]==1.0) or (i!=j and mat.values[i][j]==0.0))) return false;
}
bool check_same(matrix mat1, matrix mat2)
{
    if(mat1.row!=mat2.row or mat1.col!=mat2.col) return false;
    for(int i=0; i<mat1.row; i++)
    {
        for(int j=0; j<mat1.col; j++)
        {
            if(mat1.values[i][j]!=mat2.values[i][j])
            {
                printf("\n%f != %f\n", mat1.values[i][j], mat2.values[i][j]);
                return false;
            } 
        }
    }
    return true;
}
bool check_diagonal(matrix mat, float val)
{
    for(int i=0; i<mat.row; i++)
        for(int j=0; j<mat.col; j++)
            if(i==j and mat.values[i][j]!=val) return false;
    return true;
}
// genterater functions
struct matrix identity_matrix(int n)
{
    matrix resullt;
    resullt = onevalues(n,n,0.0);
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            if(i==j)
            {
                resullt.values[i][j] = 1;
            }
        }
    }
    return resullt;
}
struct matrix random_matrix(int row, int col, int min, int max)
{
    matrix resullt;
    resullt = onevalues(row,col,0.0);
    for(int i=0; i<row; i++)
    {
        for(int j=0; j<col; j++)
        {
            resullt.values[i][j] = randint(min, max);
        }
    }
    return resullt;
}
struct matrix augmented_matrix(matrix mat1, matrix mat2)
{
    matrix result;
    if(mat1.row!=mat2.row)
    {
        printf("Error : shape does not match");
        return result;
    }
    result = onevalues(mat1.row, mat1.col + mat2.col, 0.0);
    for(int i=0; i<result.row; i++)
    {
        for(int j=0; j<result.col; j++)
        {
            if(j<mat1.col)
            {
                result.values[i][j] = mat1.values[i][j];
            }
            else
            {
                result.values[i][j] = mat2.values[i][j-mat1.col];
            }
            
        }
    }
    return result;
}

// basic row operation (default : replace)
void basic_changerow(matrix mat, int i, int j)
{
    float temp=0;
    for(int h=0; h<mat.col; h++)
    {
        temp = mat.values[i][h];
        mat.values[i][h] = mat.values[j][h];
        mat.values[j][h] = temp;
    }
}
void basic_mulconstant(matrix mat, int i, float k)
{
    float temp=0;
    for(int h=0; h<mat.col; h++)
    {
        mat.values[i][h] = k * mat.values[i][h];
    }
}
void basic_mulandplus(matrix mat, int i, int j, float k)
{
    float temp=0;
    for(int h=0; h<mat.col; h++)
    {
        temp = k*mat.values[i][h];
        mat.values[j][h] = temp + mat.values[j][h];
    }
}
// matrix operation
struct matrix mat_constmul(matrix matt, float k)
{
    matrix mat;
    mat = copy_matrix(matt);
    for(int i=0; i<mat.row; i++)
    {
        for(int j=0; j<mat.col; j++)
        {
            mat.values[i][j] *= k;
        }
    }
    return mat;
}
struct matrix mat_sum(matrix mat1, matrix mat2)
{
    matrix result;
    result = onevalues(mat1.row, mat1.col, 0.0);
    assert(mat1.row==mat2.row and mat1.col==mat2.col);
    for(int i=0; i<mat1.row; i++)
    {
        for(int j=0; j<mat1.col; j++)
        {
            result.values[i][j] = mat1.values[i][j] + mat2.values[i][j];
        }
    }
    return result;
}
struct matrix mat_sub(matrix mat1, matrix mat2)
{
    matrix result;
    result = onevalues(mat1.row, mat1.col, 0.0);
    assert(mat1.row==mat2.row and mat1.col==mat2.col);
    for(int i=0; i<mat1.row; i++)
    {
        for(int j=0; j<mat1.col; j++)
        {
            result.values[i][j] = mat1.values[i][j] - mat2.values[i][j];
        }
    }
    return result;
}
struct matrix mat_mul(matrix mat1, matrix mat2)
{
    matrix result;
    float temp = 0.0;
    result = onevalues(mat1.row, mat2.col, 0.0);
    assert(mat1.col==mat2.row);
    for(int i=0; i<mat1.row; i++)
    {
        for(int j=0; j<mat2.col; j++)
        {
            temp = 0.0;
            for(int h=0; h<mat1.col; h++)
            {
                temp += mat1.values[i][h] * mat2.values[h][j];
            }
            result.values[i][j] = temp;
        }
    }
    return result;
}
struct matrix mat_trans(matrix mat)
{
    matrix temp;
    temp = onevalues(mat.col, mat.row, 0.0);
    for(int i=0; i<mat.row; i++)
    {
        for(int j=0; j<mat.col; j++)
        {
            temp.values[j][i] = mat.values[i][j];
        }
    }
    return temp;
}
//submatrix and cofactor
struct matrix submatrix(matrix mat, int i, int j)
{
    int ii=0;
    int jj=0;
    matrix result;
    result = onevalues(mat.row-1, mat.col-1, 0.0);
    for(int h=0; h<mat.row; h++)
    {
        if(h==i)
        {
            continue;
        }
        for(int k=0; k<mat.col; k++)
        {
            if(k==j){
                continue;
            }
            result.values[ii][jj] = mat.values[h][k];
            jj++;
        }
        ii++;
        jj=0;
    }
    return result;
}
float cofactor_recurrent(matrix mat, int i, bool dofree)
{
    if(mat.row==2)
    {
        return mat.values[0][0]*mat.values[1][1]-mat.values[1][0]*mat.values[0][1];
    }
    float result = 0.0;
    float cofac;
    for(int j=0; j<mat.col; j++)
    {
        if(mat.values[i][j]!=0)
        {
            
            cofac = n_mul(-1.0, i+j) * cofactor_recurrent(submatrix(mat, i, j), i, true);
            result += mat.values[i][j] * cofac;
        }
        
    }
    if(dofree)
    {
        free_matrix(mat);
    }
    return result;
}
float cofactor_expansion(matrix mat, int i)
{
    if(mat.row!=mat.col)
    {
        assert(false);
    }
    else if(check_zerorow(mat) or check_zerorow(mat_trans(mat)))
    {
        return 0.0;
    }
    float cofactor=0.0;
    cofactor = cofactor_recurrent(mat, i,false);
    return cofactor;
}
struct matrix cofactor_matrix(matrix mat)
{
    matrix result = onevalues(mat.row, mat.col, 0.0);
    if(mat.row!=mat.col)
    {
        assert(false);
    }
    else
    {
        for(int i=0; i<mat.row; i++)
        {
            for(int j=0; j<mat.row; j++)
            {
                result.values[i][j] = cofactor_expansion(submatrix(mat, i, j), 0) * n_mul(-1.0, i+j+2);
            }
        }
    }
    return result;
}
struct matrix adjoint_matrix(matrix mat)
{
    matrix result = cofactor_matrix(mat);
    matrix result2;
    result2 = mat_trans(result);
    free_matrix(result);
    return result2;
}
// cramer's rule
struct matrix crmsrul_answer(matrix coeffi, matrix constant)
{

    matrix answers;
    answers = onevalues(coeffi.row, 1, 0.0);
    assert(coeffi.row==coeffi.col and constant.col==1 and constant.row==coeffi.row);
    float mother, son;
    son = cofactor_expansion(coeffi, 0);
    assert(son!=0.0);
    for(int j=0; j<coeffi.col; j++)
    {
        mother = cofactor_expansion(replace_column(coeffi, constant, j), 0);
        answers.values[j][0] = mother/son;
    }
    return answers;
}
struct matrix crmsrul_inverse(matrix mat)
{
    matrix adj = adjoint_matrix(mat);
    float cofactor = cofactor_expansion(mat, 0);
    assert(cofactor!=0.0);
    adj = mat_constmul(adj, 1.0/cofactor);
    return adj;
}
// elimination
bool is_REF(matrix mat)
{
    int temp=-1;
    bool zerorow;
    bool bzerorow=false;
    for(int i=0; i<mat.row; i++)
    {
        zerorow=true;
        for(int j=0; j<mat.col; j++)
        {
            if(mat.values[i][j]!=0.0)
            {
                zerorow=false;
                if(mat.values[i][j]!=1.0)
                {
                    return false;
                }
                else
                {
                    if(temp < j)
                    {
                        temp=j;
                        break;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        if(!zerorow and bzerorow)
        {
            return false;
        }
        bzerorow = zerorow;
    }
    return true;
}
struct matrix REF(matrix mat)
{
    int first, ii;
    float constt;
    matrix result = copy_matrix(mat);
    for(int j=0; j<result.col; j++)
    {
        first=-1;
        ii=j;
        while(ii<result.row)
        {
            if(result.values[ii][j]!=0.0)
            {
                first=ii;
                break;
            }
            ii++;
        }
        if(first==-1) continue;
        basic_mulconstant(result, first, 1/result.values[first][j]);
        for(int i=j; i<result.row; i++)
        {
            if(i==first) continue;
            constt = -result.values[i][j]/result.values[first][j];
            basic_mulandplus(result, first, i, constt);
        }
        basic_changerow(result, first, j);
    }
    return result;
}
struct matrix RREF(matrix mat)
{
    int first, ii;
    float constt;
    matrix result = copy_matrix(mat);
    for(int j=0; j<result.col; j++)
    {
        first=-1;
        ii=j;
        while(ii<result.row)
        {
            if(result.values[ii][j]!=0.0)
            {
                first=ii;
                break;
            }
            ii++;
        }
        if(first==-1) continue;
        basic_mulconstant(result, first, 1/result.values[first][j]);
        for(int i=0; i<result.row; i++)
        {
            if(i==first) continue;
            constt = -result.values[i][j]/result.values[first][j];
            basic_mulandplus(result, first, i, constt);
        }
        basic_changerow(result, first, j);
    }
    return result;
}
struct matrix elimin_inverse(matrix mat)
{
    assert(mat.row==mat.col);
    matrix aug = augmented_matrix(mat, identity_matrix(mat.row));
    matrix result = RREF(aug);
    matrix iden = clip_by_col(result, 0, mat.row);
    matrix inverse = clip_by_col(result, mat.row, mat.row*2);
    free_matrix(aug);
    free_matrix(result);
    if(check_identity(tidy_matrix(iden)))
    {
        free_matrix(iden);
        return inverse;
    }
    else
    {
        assert(false);
    }
}
struct matrix RREF_answer(matrix coeffi, matrix constant)
{
    assert(coeffi.row==coeffi.col and constant.col==1 and constant.row==coeffi.row);
    matrix answers;
    matrix new_coeffi;
    
    matrix aug = augmented_matrix(coeffi, constant);
    matrix rref = RREF(aug);
    answers = clip_by_col(rref, coeffi.col, coeffi.col+1);
    new_coeffi = clip_by_col(rref, 0, coeffi.col);
    assert(check_identity(tidy_matrix(new_coeffi)));

    free_matrix(aug);
    free_matrix(rref);
    return answers;
}
struct matrix REF_answer(matrix coeffi, matrix constant)
{
    assert(coeffi.row==coeffi.col and constant.col==1 and constant.row==coeffi.row);
    matrix answers;
    matrix aug = augmented_matrix(coeffi, constant);
    matrix rref = REF(aug);
    matrix new_coeffi = clip_by_col(rref, 0, coeffi.col);
    matrix new_constant = clip_by_col(rref, coeffi.col, coeffi.col+1);

    int valnum = 0;
    float mem=0.0;
    if(check_diagonal(tidy_matrix(new_coeffi), 1.0) and !check_zerorow(tidy_matrix(new_coeffi)))
    {
        // i는 행 
        for(int i=coeffi.row-2; i>-1; i--)
        {
            //.valnum 변수의 개수 
            valnum = coeffi.row - i;
            mem=0.0;
            for(int j=coeffi.col-valnum+1; j<coeffi.col; j++)
            {
                mem += new_coeffi.values[i][j]*new_constant.values[j][0];
                new_coeffi.values[i][j] = 0;
            }
            new_constant.values[i][0] -= mem;
        }
    }
    else
    {
        printf("REF_answer : Error : no answer");
        assert(false);
    }

    free_matrix(aug);
    free_matrix(rref);
    free_matrix(new_coeffi);
    return new_constant;
}


int main()
{
    srand(time(NULL));

}