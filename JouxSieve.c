#define STB_DS_IMPLEMENTATION
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpzi.h>
#include <flint/fmpq.h>
#include <flint/fmpz_factor.h>
#include "stb_ds.h"
#define LINE_BUF_SIZE 4096
#define INFINITE_LOOP_CHECK_SUM 100000

//clear && gcc JouxSieve.c -o m.o -lm -lgmp -lmpfr -lflint && ./m.o 
typedef struct hash_table_struct *HashTable;
typedef struct hash_table_entry_struct *HashTableEntry;
struct hash_table_entry_struct
{
	fmpz_t key;
	fmpz_t value;
};
struct hash_table_struct
{
	HashTableEntry **entry;
	int *primeList;
	fmpz_t ell;
	fmpz_t modulo;
	fmpz_t unknownBase;
	fmpz_t unknownBaseProjection;
	fmpz_t generatorPmin1;
	size_t capacity;
};
HashTableEntry CreateEntry();
void FreeHashTableEntry(HashTableEntry entry);
HashTable CreateHashTable(size_t capacity, char *unknownBaseString, char *unknownBaseProjectionString, char *ellString, char *generatorPmin1String, int base);
bool FindElementHashTable(HashTable table, fmpz_t search, fmpz_t temp0, fmpz_t result);
void DestroyHashTable(HashTable table);
HashTableEntry CreateEntry()
{
	HashTableEntry entry = malloc(sizeof(struct hash_table_entry_struct));
	fmpz_init(entry->key);
	fmpz_init(entry->value);
	return entry;
}
void FreeHashTableEntry(HashTableEntry entry)
{
	if(entry)
	{
		fmpz_clear(entry->key);
		fmpz_clear(entry->value);
		free(entry);
	}
}
HashTable CreateHashTable(size_t capacity, char *unknownBaseString, char *unknownBaseProjectionString, char *ellString, char *generatorPmin1String, int base)
{
	HashTable table = malloc(sizeof(struct hash_table_struct));
	table->primeList = NULL;
	fmpz_init(table->modulo);fmpz_init(table->unknownBase);fmpz_init(table->unknownBaseProjection);fmpz_init(table->ell);
	fmpz_init(table->generatorPmin1);
	assert(capacity > 0);
	fmpz_set_ui(table->modulo, capacity);
	fmpz_set_str(table->unknownBase, unknownBaseString, base);
	fmpz_set_str(table->unknownBaseProjection, unknownBaseProjectionString, base);
	fmpz_set_str(table->ell, ellString, base);
	fmpz_set_str(table->generatorPmin1, generatorPmin1String, base);
	table->capacity = capacity;
	table->entry = malloc(capacity * sizeof(HashTableEntry*));
	for(size_t i = 0; i < capacity; i++)
	{
		table->entry[i] = NULL;
	}
	return table;
}
void DestroyHashTable(HashTable table)
{
	for(size_t i = 0; i < table->capacity; i++)
	{
		for(size_t j = 0; j < arrlen(table->entry[i]); j++)
		{
			FreeHashTableEntry(table->entry[i][j]);
		}
		arrfree(table->entry[i]);
	}
	fmpz_clear(table->generatorPmin1);
	fmpz_clear(table->ell);
	fmpz_clear(table->unknownBaseProjection);
	fmpz_clear(table->unknownBase);
	fmpz_clear(table->modulo);
	arrfree(table->primeList);
	free(table->entry);
	free(table);
}
void PrintTableContents(HashTable table)
{
	for(size_t i = 0; i < table->capacity; i++)
	{
		if(arrlen(table->entry[i]) > 0)
		{
			printf("Mod %3ld\n", i);
			for(size_t j = 0; j < arrlen(table->entry[i]); j++)
			{
				printf("\t");
				fmpz_print(table->entry[i][j]->key);
				printf(" : ");
				fmpz_print(table->entry[i][j]->value);
				printf("\n");
			}
		}
	}
}
bool FindElementHashTable(HashTable table, fmpz_t search, fmpz_t temp0, fmpz_t result)
{
	bool found = false;
	//Find modulo
	fmpz_mod(temp0, search, table->modulo);
	int searchIndex = fmpz_get_ui(temp0);
	for(size_t i = 0; i < arrlen(table->entry[searchIndex]); i++)
	{
		if(fmpz_cmp(search, table->entry[searchIndex][i]->key) == 0)
		{
			fmpz_set(result, table->entry[searchIndex][i]->value);
			found = true;
			break;
		}
	}
	return found;
}

void SetElementHashTable(HashTable table, fmpz_t n, fmpz_t temp0, fmpz_t result)
{
	//Find modulo
	fmpz_mod(temp0, n, table->modulo);
	int searchIndex = fmpz_get_ui(temp0);
	HashTableEntry entry = CreateEntry();
	fmpz_set(entry->key, n);
	fmpz_set(entry->value, result);
	arrput(table->entry[searchIndex],entry);
}

void Positive_GaussReduce(fmpz_t x1, fmpz_t y1, fmpz_t x2, fmpz_t y2)
{
	fmpz_t mu, dot12, dot11,temp0, tmp2, tmp3;
	fmpz_init(mu);fmpz_init(dot12);fmpz_init(dot11);fmpz_init(temp0);fmpz_init(tmp2);fmpz_init(tmp3);
	int infiniteLoopCheck = 0;
	while(1)
	{
		/*dot12 = v1·v2 */
		fmpz_mul(dot12, x1, x2);fmpz_mul(temp0, y1, y2);fmpz_add(dot12, dot12, temp0);
		/*dot11 = v1·v1 */
		fmpz_mul(dot11, x1, x1);fmpz_mul(temp0, y1, y1);fmpz_add(dot11, dot11, temp0);
		/*avoid division by zero */
		if(fmpz_is_zero(dot11)){break;}
		/*mu = floor((v1·v2)/(v1·v1))*/
		fmpz_fdiv_q(mu, dot12, dot11);
		/*v2 = v2 - mu*v1 */
		fmpz_mul(temp0, mu, x1);fmpz_sub(x2, x2, temp0);fmpz_mul(temp0, mu, y1);fmpz_sub(y2, y2, temp0);
		/*||v1||^2 */
		fmpz_mul(temp0, x1, x1);fmpz_mul(tmp2, y1, y1);fmpz_add(temp0, temp0, tmp2);
		/*||v2||^2 */
		fmpz_mul(tmp2, x2, x2);fmpz_mul(tmp3, y2, y2);fmpz_add(tmp2, tmp2, tmp3);
		/*If ||v2|| < ||v1|| then swap*/
		if(fmpz_cmp(tmp2, temp0) < 0){fmpz_swap(x1, x2);fmpz_swap(y1, y2);}
		else{break;}
		infiniteLoopCheck++;
		if(infiniteLoopCheck > INFINITE_LOOP_CHECK_SUM){printf("Gauss reduction may be in infinite loop\n");abort();}
	}
	fmpz_clear(mu);fmpz_clear(dot12);fmpz_clear(dot11);fmpz_clear(temp0);fmpz_clear(tmp2);fmpz_clear(tmp3);
}

void PrintBasis2D(fmpz_t x1, fmpz_t y1, fmpz_t x2, fmpz_t y2)
{
	fmpz_print(x1);printf(" ");fmpz_print(x2);printf("\n");
	fmpz_print(y1);printf(" ");fmpz_print(y2);printf("\n");
}

bool ProcessLogFileLine(char *line, fmpz_t primeHolder, fmpz_t logHolder)
{
	bool goodLine = true;
	
	char *tokens[8];
	int count = 0;
	char *tok = strtok(line, " \t\r\n");
	while(tok && count < 8){tokens[count++] = tok;tok = strtok(NULL, " \t\r\n");}
	if(count < 4){goodLine = false;return goodLine;}
	char *index = tokens[0];
	char *p = tokens[1];
	char *side = tokens[2];
	char *r = tokens[3];
	if(strcmp(p, "bad") == 0 && strcmp(side, "ideals") == 0){goodLine = false;return goodLine;}
	if(strcmp(p, "SM") == 0){goodLine = false;return goodLine;}
	int sideValue = atoi(side);
	if(strcmp(r, "rat") == 0 && sideValue == 0)
	{
		char *logString = tokens[4];
		int64_t primeValue = strtoll(p, NULL, 16);
		fmpz_set_str(primeHolder, p, 16);
		fmpz_set_str(logHolder,logString, 10);
		//printf("%lu %s\n",primeValue, logString);
		return goodLine;
	}
	goodLine = false;return goodLine;
}
void LoadLogFile(HashTable logTable,char *logFileName)
{
	FILE *fp = fopen(logFileName, "r");assert(fp != NULL);
	char line[LINE_BUF_SIZE];
	fmpz_t primeHolder, logHolder,temp0,temp1,targetBaseInverse;
	fmpz_init(primeHolder);fmpz_init(logHolder);fmpz_init(temp0);fmpz_init(temp1);fmpz_init(targetBaseInverse);
	fmpz_invmod(targetBaseInverse, logTable->unknownBaseProjection,logTable->ell);
	//Read first line
	fgets(line, sizeof(line), fp);
	//Read other lines
	while(fgets(line, sizeof(line), fp))
	{
		bool validLine = ProcessLogFileLine(line, primeHolder, logHolder);
		if(validLine)
		{
			//printf("p:");fmpz_print(primeHolder);printf(", ");printf("l:");fmpz_print(logHolder);printf("\n");
			bool found = FindElementHashTable(logTable, primeHolder, temp0, temp1);
			int currentPrime = fmpz_get_ui(primeHolder);
			arrput(logTable->primeList, currentPrime);
			assert(found == false);
			SetElementHashTable(logTable, primeHolder, temp0, logHolder);	
		}
		//break;
	}
	fclose(fp);
	fmpz_clear(primeHolder);fmpz_clear(logHolder);fmpz_clear(temp0);fmpz_clear(temp1);fmpz_clear(targetBaseInverse);
}

/*Test target = A/B mod p*/
bool TestRationalReconstruction(fmpz_t target, fmpz_t prime, fmpz_t A, fmpz_t B)
{
	bool validReconstruction = false;
	fmpz_t temp0;fmpz_init(temp0);
	fmpz_invmod(temp0,B,prime);
	fmpz_mul(temp0, temp0, A);
	fmpz_mod(temp0,temp0,prime);
	//printf("A:");fmpz_print(A);printf("\n");
	//printf("B:");fmpz_print(B);printf("\n");
	//printf("temp0:");fmpz_print(temp0);printf("\n");
	validReconstruction = fmpz_cmp(temp0, target) == 0;
	fmpz_clear(temp0);
	return validReconstruction;
}

/*Test identity target = (k1A1 + k2A2) / (k1B1 + k2B2) mod p*/
bool TestRationalIdentity(fmpz_t target, fmpz_t prime, fmpz_t A1, fmpz_t B1, fmpz_t A2, fmpz_t B2,flint_rand_t rngState)
{
	bool validReconstruction = true;
	fmpz_t temp0,num,den,k1,k2;
	fmpz_init(temp0);fmpz_init(num);fmpz_init(den);fmpz_init(k1);fmpz_init(k2);
	int randomTestCount = 500;
	for(int i = 0; i < randomTestCount; i++)
	{
		fmpz_randm(k1, rngState, prime);fmpz_randm(k2, rngState, prime);
		fmpz_mul(temp0,k1,A1);fmpz_mul(num,k2,A2);fmpz_add(num,num,temp0);fmpz_mod(num, num, prime);
		fmpz_mul(temp0,k1,B1);fmpz_mul(den,k2,B2);fmpz_add(den,den,temp0);fmpz_mod(den, den, prime);
		fmpz_invmod(temp0, den, prime);fmpz_mul(temp0, temp0, num);fmpz_mod(temp0, temp0, prime);
		//printf("Z:");fmpz_print(temp0);printf("\n");
		if(fmpz_cmp(temp0, target) != 0)
		{
			validReconstruction = false;
			break;
		}
	}
	
	fmpz_clear(temp0);fmpz_clear(num);fmpz_clear(den);fmpz_clear(k1);fmpz_clear(k2);
	return validReconstruction;
}
void TestLoadLogs()
{
	char *logFileName = "../Data/p95.dlog";
	size_t capacity = 100003;
	char *unknownBaseString           = "101259726400208184076718506974609171605435438168191661079316293670904622880200";
	char *unknownBaseProjectionString = "97890684989265514181240123615056668210036264581619234259585154409438989672892";
	char *ellString        = "205115282021455665897114700593932402728804164701536103180137503955397371";
	char *generatorPmin1String = "3";
	int base = 10;
	HashTable logTable = CreateHashTable(capacity,unknownBaseString,unknownBaseProjectionString, ellString, generatorPmin1String,base);
	LoadLogFile(logTable, logFileName);
	PrintTableContents(logTable);
	DestroyHashTable(logTable);
}

void TestSieveIdentity()
{
	char *logFileName = "../Data/p95.dlog";
	size_t capacity = 100003;
	char *unknownBaseString           = "101259726400208184076718506974609171605435438168191661079316293670904622880200";
	char *unknownBaseProjectionString = "97890684989265514181240123615056668210036264581619234259585154409438989672892";
	char *ellString                   = "205115282021455665897114700593932402728804164701536103180137503955397371";
	char *primeString                 = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
	char *generatorPmin1String = "3";int base = 10;
	HashTable logTable = CreateHashTable(capacity,unknownBaseString,unknownBaseProjectionString, ellString, generatorPmin1String,base);
	LoadLogFile(logTable, logFileName);
	
	int testCount = 10000;
	flint_rand_t rngState;flint_rand_init(rngState);
	fmpz_t currentPower,prime,x1, y1, x2, y2;
	fmpz_init(currentPower);fmpz_init(prime);fmpz_init(x1);fmpz_init(y1);fmpz_init(x2);fmpz_init(y2);
	fmpz_set(currentPower, logTable->generatorPmin1);fmpz_set_str(prime, primeString, 16);
	for(int i = 1; i < testCount; i++)
	{
		fmpz_set(x1, currentPower);fmpz_set_ui(y1, 1);
		fmpz_set(x2, prime);fmpz_set_ui(y2, 0);
		Positive_GaussReduce(x1, y1, x2, y2);
		bool valid1 = TestRationalReconstruction(currentPower, prime, x1, y1);assert(valid1 == true);
		bool valid2 = TestRationalReconstruction(currentPower, prime, x2, y2);assert(valid2 == true);
		bool valid3 = TestRationalIdentity(currentPower, prime, x1, y1, x2, y2,rngState);assert(valid3 == true);
		printf("\n3^(%d):", i);fmpz_print(currentPower);printf("\n");
		PrintBasis2D(x1, y1, x2, y2);
		fmpz_mul(currentPower, currentPower, logTable->generatorPmin1);
		fmpz_mod(currentPower,currentPower, prime);
	}
	flint_rand_clear(rngState);
	DestroyHashTable(logTable);
	fmpz_clear(currentPower);fmpz_clear(prime);fmpz_clear(x1);fmpz_clear(y1);fmpz_clear(x2);fmpz_clear(y2);
}

int SinglePassSieve(int *largestK1Holder,int k2, size_t k1Max, int *k1Array,int *primeList,fmpz_t a, fmpz_t c, fmpz_t b, fmpz_t d)
{
	int largestK1 = 0;
	int largestK1Index = 0;
	fmpz_t temp0,currentPrime;
	fmpz_init(temp0);fmpz_init(currentPrime);
	for(size_t primeIndex = 0; primeIndex < arrlen(primeList); primeIndex++)
	{
		fmpz_set_ui(currentPrime, primeList[primeIndex]);
		//Test gcd
		if(fmpz_invmod(temp0, a, currentPrime) != 0)
		{
			//Find first valid index
			fmpz_mul(temp0, temp0, b);
			fmpz_mul_si(temp0, temp0, -k2);
			fmpz_mod(temp0, temp0, currentPrime);
			size_t validIndex = fmpz_get_ui(temp0);
			for(size_t currentKIndex = validIndex; currentKIndex < k1Max; currentKIndex += primeList[primeIndex])
			{
				k1Array[currentKIndex] += (int)(1000 * log(primeList[primeIndex]));
				if(k1Array[currentKIndex] > largestK1)
				{
					largestK1 = k1Array[currentKIndex];
					largestK1Index = currentKIndex;
				}
			}
			//printf("Valid0: %ld mod %d\n",validIndex,primeList[primeIndex]);
		}
		if(fmpz_invmod(temp0, c, currentPrime) != 0)
		{
			//Find first valid index
			fmpz_mul(temp0, temp0, d);
			fmpz_mul_si(temp0, temp0, -k2);
			fmpz_mod(temp0, temp0, currentPrime);
			size_t validIndex = fmpz_get_ui(temp0);
			for(size_t currentKIndex = validIndex; currentKIndex < k1Max; currentKIndex += primeList[primeIndex])
			{
				k1Array[currentKIndex] += (int)(1000 * log(primeList[primeIndex]));
				if(k1Array[currentKIndex] > largestK1)
				{
					largestK1 = k1Array[currentKIndex];
					largestK1Index = currentKIndex;
				}
			}
			//printf("Valid1: %ld mod %d\n",validIndex,primeList[primeIndex]);
		}
	}	
	fmpz_clear(temp0);fmpz_clear(currentPrime);
	*largestK1Holder = largestK1;
	return largestK1Index;
}

void PrintFactors(fmpz_factor_t factors)
{
	for(slong i = 0; i < factors->num; i++)
	{
		int size = fmpz_sizeinbase(&factors->p[i], 2);
		printf("(|%d| ", size);
		fmpz_print(&factors->p[i]); printf(" ^ ");
		fmpz_print(&factors->exp[i]); printf("), ");
	}
	printf("\n");
}

bool TestFactorization(fmpz_t target, fmpz_t prime, fmpz_t A1, fmpz_t B1, fmpz_t A2, fmpz_t B2, int k1, int k2)
{
	bool validReconstruction = true;
	fmpz_t temp0,num,den;
	fmpz_init(temp0);fmpz_init(num);fmpz_init(den);
	int randomTestCount = 500;
	
	fmpz_mul_ui(temp0,A1,k1);fmpz_mul_ui(num,A2,k2);fmpz_add(num,num,temp0);fmpz_mod(num, num, prime);
	fmpz_mul_ui(temp0,B1,k1);fmpz_mul_ui(den,B2,k2);fmpz_add(den,den,temp0);//fmpz_mod(den, den, prime);
	fmpz_invmod(temp0, den, prime);fmpz_mul(temp0, temp0, num);fmpz_mod(temp0, temp0, prime);
	if(fmpz_cmp(temp0, target) != 0)
	{
		validReconstruction = false;
	}
	fmpz_factor_t factorNum, factorDen;
	fmpz_factor_init(factorNum);fmpz_factor_init(factorDen);
	
	int proved = 1;
	int factorNumResult = fmpz_factor_smooth(factorNum, num, 20, proved);
	int factorDenResult = fmpz_factor_smooth(factorDen, den, 20, proved);
	printf("Num:");fmpz_print(num);printf("\n");PrintFactors(factorNum);
	printf("Den:");fmpz_print(den);PrintFactors(factorDen);printf("\n\n");
	fmpz_factor_clear(factorNum);fmpz_factor_clear(factorDen);
	fmpz_clear(temp0);fmpz_clear(num);fmpz_clear(den);
	return validReconstruction;
}

bool Sieve_k1_k2(fmpz_t target, fmpz_t prime, fmpz_t A1, fmpz_t B1, fmpz_t A2, fmpz_t B2, int *primeList)
{
	bool foundGoodPair = false;
	size_t k1Max = 2000000;
	int k2Max = 16044;
	int largestK1Index_pre = 0;
	//Fix different k2's
	int *k1Array = calloc(k1Max, sizeof(int));
	for(int k2 = 16043; k2 < k2Max; k2++)
	{
		if(k2 == 0){continue;}
		//Reset k1's
		memset(k1Array, 0, k1Max * sizeof(int));
		int largestK1Holder = 0;
		int largestK1Index = SinglePassSieve(&largestK1Holder, k2, k1Max, k1Array, primeList, A1, B1, A2, B2);
		printf("Found pair: (%7d, %7d):%5d %5d \n",largestK1Index, k2,largestK1Index_pre,largestK1Holder);
		if(largestK1Holder > largestK1Index_pre)
		{
			bool valid = TestFactorization(target, prime, A1, B1, A2, B2, largestK1Index, k2);
			largestK1Index_pre = largestK1Holder;
		}
		
	}
	free(k1Array);
	return foundGoodPair;
}

void TestSieveArray()
{
	char *logFileName = "Data/p95.dlog";
	size_t capacity = 100003;
	char *unknownBaseString           = "101259726400208184076718506974609171605435438168191661079316293670904622880200";
	char *unknownBaseProjectionString = "97890684989265514181240123615056668210036264581619234259585154409438989672892";
	char *ellString                   = "205115282021455665897114700593932402728804164701536103180137503955397371";
	char *primeString                 = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
	char *generatorPmin1String = "3";int base = 10;
	HashTable logTable = CreateHashTable(capacity,unknownBaseString,unknownBaseProjectionString, ellString, generatorPmin1String,base);
	LoadLogFile(logTable, logFileName);
	
	int testCount = 2;
	int powStart = 10000;
	flint_rand_t rngState;flint_rand_init(rngState);
	fmpz_t currentPower,prime,x1, y1, x2, y2;
	fmpz_init(currentPower);fmpz_init(prime);fmpz_init(x1);fmpz_init(y1);fmpz_init(x2);fmpz_init(y2);
	fmpz_set_str(prime, primeString, 16);
	fmpz_powm_ui(currentPower, logTable->generatorPmin1, powStart,prime);
	
	for(int i = 1; i < testCount; i++)
	{
		fmpz_set(x1, currentPower);fmpz_set_ui(y1, 1);
		fmpz_set(x2, prime);fmpz_set_ui(y2, 0);
		Positive_GaussReduce(x1, y1, x2, y2);
		printf("\n3^(%d):", i);fmpz_print(currentPower);printf("\n");
		PrintBasis2D(x1, y1, x2, y2);
		bool foundGoodk1k2  = Sieve_k1_k2(currentPower, prime, x1, y1, x2, y2, logTable->primeList);
		
		
		fmpz_mul(currentPower, currentPower, logTable->generatorPmin1);
		fmpz_mod(currentPower,currentPower, prime);
	}
	flint_rand_clear(rngState);
	DestroyHashTable(logTable);
	fmpz_clear(currentPower);fmpz_clear(prime);fmpz_clear(x1);fmpz_clear(y1);fmpz_clear(x2);fmpz_clear(y2);
}


int main()
{
	TestSieveArray();
	//TestLoadLogs();
	//TestSieveIdentity();
	flint_cleanup();
	return 0;
}
