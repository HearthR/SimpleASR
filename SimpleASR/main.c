#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmm.h"

wordHmmType word_hmm[N_WORD];

void buildWordHMM()
{
	FILE *dict_fp;
	errno_t err = fopen_s(&dict_fp, "dictionary.txt", "rt");
	if (err == 0)
	{
		printf("File open successed!\n");
	}
	else
	{
		printf("File open failed!\n");
		return;
	}

	char buffer[25];
	char buffer_t[25];
	int i = 0;
	char* token;
	int dim = 0;
	
	while (fgets(buffer, sizeof(buffer), dict_fp) != NULL)
	{
		strcpy(buffer_t, buffer);

		token = strtok(buffer_t, " \t");
		strcpy(&word_hmm[i].name, token);
		printf("%s : ", word_hmm[i].name);

		token = strtok(NULL, " \n");
		while (token != NULL)
		{
			if (strcmp(token, "sp") == 0)
			{
				dim++;
			}
			else
			{
				dim += 3;
			}
			token = strtok(NULL, " \n");
			// tp를 어떻게 만들지. -> ey t sp일때 tp 차원은 1 + 3 + 3 + 1 + 1 로 hmm.h처럼 만들까?ㄴ
		}
		printf("%d\n", dim);
		word_hmm[i].tp_dim = dim + 2;

		word_hmm[i].tp = calloc(1, (sizeof(float) * (dim + 2)) * (sizeof(float) * (dim + 2)));
		
		dim = 0;
		hmmType* hmm_temp;
		token = strtok(buffer, " \t");
		token = strtok(NULL, " \n");
		while (token != NULL)
		{
			for (int x = 0; x < 21; x++)
			{
				if (strcmp(phones[x].name, token) == 0)
				{
					hmm_temp = &phones[x];
					break;
				}
			}

			if (strcmp(token, "sil") == 0)
			{
				for (int x = 0; x < 5; x++)
				{
					for (int y = 0; y < 5; y++)
					{
						word_hmm[i].tp[x * 5 + y] = hmm_temp->tp[x][y];
					}
				}
			}
			else if (strcmp(token, "sp") == 0)
			{
				word_hmm[i].tp[word_hmm[i].tp_dim * dim + dim + 2] = word_hmm[i].tp[word_hmm[i].tp_dim * dim + dim + 1] * hmm_temp->tp[0][2];
				word_hmm[i].tp[word_hmm[i].tp_dim * dim + dim + 1] *= hmm_temp->tp[0][1];
				word_hmm[i].tp[word_hmm[i].tp_dim * (dim + 1) + dim + 1] = hmm_temp->tp[1][1];
				word_hmm[i].tp[word_hmm[i].tp_dim * (dim + 1) + dim + 2] = hmm_temp->tp[1][2];
			}
			else
			{
				if (dim == 0)
				{
					for (int x = 0; x < 4; x++)
					{
						word_hmm[i].tp[word_hmm[i].tp_dim * (dim + x) + dim + x + 1] = hmm_temp->tp[x][x + 1];
						word_hmm[i].tp[word_hmm[i].tp_dim * (dim + x + 1) + dim + x + 1] = hmm_temp->tp[x + 1][x + 1];
					}
				}
				else
				{
					word_hmm[i].tp[word_hmm[i].tp_dim * (dim) + dim + 1] *= hmm_temp->tp[0][1];
					word_hmm[i].tp[word_hmm[i].tp_dim * (dim + 1) + dim + 1] = hmm_temp->tp[1][1];
					for (int x = 1; x < 4; x++)
					{
						word_hmm[i].tp[word_hmm[i].tp_dim * (dim + x) + dim + x + 1] = hmm_temp->tp[x][x + 1];
						word_hmm[i].tp[word_hmm[i].tp_dim * (dim + x + 1) + dim + x + 1] = hmm_temp->tp[x + 1][x + 1];
					}
				}
				dim += 3;
			}
			
			printf("%s\n", token);
			token = strtok(NULL, " \n");
		}

		for (int x = 0; x < word_hmm[i].tp_dim; x++)
		{
			for (int y = 0; y < word_hmm[i].tp_dim; y++)
			{
				printf("%.5f   ", word_hmm[i].tp[word_hmm[i].tp_dim * x + y]);
			}
			printf("\n");
		}

		printf("\n\n\n");
		i++;
		dim = 0;
	}


	
}



int main()
{
	buildWordHMM();





	return 0;
}