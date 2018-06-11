#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "hmm.h"

#define _USE_MATH_DEFINES
#include <math.h>

wordHmmType word_hmm[N_WORD];
univHmmType univ_hmm;
unigramType unigram;


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
		word_hmm[i].phones_num = 0;

		token = strtok(NULL, " \n");
		while (token != NULL)
		{
			if (strcmp(token, "sp") == 0)
			{
				dim++;
				word_hmm[i].phones_num++;
			}
			else
			{
				dim += 3;
				word_hmm[i].phones_num++;
			}
			token = strtok(NULL, " \n");
		}
		word_hmm[i].tp_dim = dim + 2;

		word_hmm[i].tp = calloc(1, (sizeof(float) * (dim + 2)) * (sizeof(float) * (dim + 2)));
		
		dim = 0;
		hmmType* hmm_temp;
		int ph_idx = 0;
		token = strtok(buffer, " \t");
		token = strtok(NULL, " \n");
		while (token != NULL)
		{
			for (int x = 0; x < 21; x++)
			{
				if (strcmp(phones[x].name, token) == 0)
				{
					hmm_temp = &phones[x];
					word_hmm[i].phones_index[ph_idx] = x;
					ph_idx++;
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
			else if (strcmp(token, "zero") == 0)
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
					word_hmm[i].tp[word_hmm[i].tp_dim * (dim)+dim + 1] *= hmm_temp->tp[0][1];
					word_hmm[i].tp[word_hmm[i].tp_dim * (dim + 1) + dim + 1] = hmm_temp->tp[1][1];
					for (int x = 1; x < 4; x++)
					{
						word_hmm[i].tp[word_hmm[i].tp_dim * (dim + x) + dim + x + 1] = hmm_temp->tp[x][x + 1];
						word_hmm[i].tp[word_hmm[i].tp_dim * (dim + x + 1) + dim + x + 1] = hmm_temp->tp[x + 1][x + 1];
					}
				}
				dim += 3;
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
			
			token = strtok(NULL, " \n");
		}

		word_hmm[i].out_idx1 = word_hmm[i].tp_dim * (word_hmm[i].tp_dim - 3) + word_hmm[i].tp_dim - 1;
		word_hmm[i].out_idx2 = word_hmm[i].tp_dim * (word_hmm[i].tp_dim - 2) + word_hmm[i].tp_dim - 1;

		i++;
		dim = 0;
	}

	fclose(dict_fp);
}

void buildUniversalHMM()
{
	int univ_dim = 0;
	FILE *bigram_fp;
	errno_t err = fopen_s(&bigram_fp, "bigram.txt", "rt");
	if (err == 0)
	{
		printf("File open successed!\n");
	}
	else
	{
		printf("File open failed!\n");
		return;
	}

	for (int i = 0; i < N_WORD; i++)
	{
		univ_dim += word_hmm[i].tp_dim - 2;
	}
	univ_dim += 2;
	univ_hmm.tp_dim = univ_dim;

	univ_hmm.tp = calloc(1, (sizeof(float)*univ_hmm.tp_dim) * (sizeof(float)*univ_hmm.tp_dim));

	int acc_dim = 0;
	int tp_phone_i = 1;
	for (int i = 0; i < N_WORD; i++)
	{
		strcpy(univ_hmm.word_list[i].name, word_hmm[i].name);
		univ_hmm.word_list[i].in_index = acc_dim + 1;
		for (int j = 1; j < word_hmm[i].tp_dim - 1; j++)
		{
			for (int k = 1; k < word_hmm[i].tp_dim - 1; k++)
			{
				univ_hmm.tp[univ_hmm.tp_dim * (acc_dim + j) + acc_dim + k] = word_hmm[i].tp[word_hmm[i].tp_dim * j + k];
			}
		}
		univ_hmm.word_list[i].state_num = word_hmm[i].tp_dim - 2;
		int st_tmp = 0;
		for (int j = 0; j < word_hmm[i].phones_num; j++)
		{
			if (strcmp(phones[word_hmm[i].phones_index[j]].name, "sp") == 0)
			{
				univ_hmm.word_list[i].state_index[st_tmp][0] = word_hmm[i].phones_index[j];
				univ_hmm.word_list[i].state_index[st_tmp][1] = 0;
				univ_hmm.tp_phones[tp_phone_i] = word_hmm[i].phones_index[j];
				st_tmp++;
				tp_phone_i++;
			}
			else
			{
				for (int x = 0; x < N_STATE; x++)
				{
					univ_hmm.word_list[i].state_index[st_tmp][0] = word_hmm[i].phones_index[j];
					univ_hmm.word_list[i].state_index[st_tmp][1] = x;
					univ_hmm.tp_phones[tp_phone_i] = word_hmm[i].phones_index[j];
					st_tmp++;
					tp_phone_i++;
				}
			}
		}
		acc_dim += word_hmm[i].tp_dim - 2;
		univ_hmm.word_list[i].out_index = acc_dim;
	}

	char before[6];
	char after[6];
	float b_prob;
	int before_idx, after_idx;
	int tp_idx;
	int zero1_idx = -1, zero2_idx = -1;

	for (int i = 0; i < N_WORD; i++)
	{
		if (strcmp("zero", univ_hmm.word_list[i].name) == 0)
		{
			if (zero1_idx == -1)
			{
				zero1_idx = i;
			}
			else
			{
				zero2_idx = i;
			}
		}
	}
	
	while (!feof(bigram_fp))
	{
		fscanf(bigram_fp, "%s\t%s\t%f\n", before, after, &b_prob);

		for (int i = 0; i < N_WORD; i++)
		{
			if (strcmp(before, univ_hmm.word_list[i].name) == 0)
			{
				before_idx = i;
			}
			if (strcmp(after, univ_hmm.word_list[i].name) == 0)
			{
				after_idx = i;
			}
		}

		if (strcmp(before, "zero") == 0)
		{
			if (strcmp(after, "zero") == 0)
			{
				tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[zero1_idx].out_index - 1) + univ_hmm.word_list[zero2_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero1_idx].tp[word_hmm[zero1_idx].out_idx1] * b_prob;
				tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[zero1_idx].out_index + univ_hmm.word_list[zero2_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero1_idx].tp[word_hmm[zero1_idx].out_idx2] * b_prob;

				tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[zero2_idx].out_index - 1) + univ_hmm.word_list[zero1_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero2_idx].tp[word_hmm[zero2_idx].out_idx1] * b_prob;
				tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[zero2_idx].out_index + univ_hmm.word_list[zero1_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero2_idx].tp[word_hmm[zero2_idx].out_idx2] * b_prob;
			}
			else
			{
				tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[zero1_idx].out_index - 1) + univ_hmm.word_list[after_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero1_idx].tp[word_hmm[zero1_idx].out_idx1] * b_prob;
				tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[zero1_idx].out_index + univ_hmm.word_list[after_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero1_idx].tp[word_hmm[zero1_idx].out_idx2] * b_prob;

				tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[zero2_idx].out_index - 1) + univ_hmm.word_list[after_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero2_idx].tp[word_hmm[zero2_idx].out_idx1] * b_prob;
				tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[zero2_idx].out_index + univ_hmm.word_list[after_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero2_idx].tp[word_hmm[zero2_idx].out_idx2] * b_prob;
			}
		}
		else if (strcmp(after, "zero") == 0)
		{
			if (strcmp(before, "zero") == 0)
			{
				tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[zero1_idx].out_index - 1) + univ_hmm.word_list[zero2_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero1_idx].tp[word_hmm[zero1_idx].out_idx1] * b_prob;
				tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[zero1_idx].out_index + univ_hmm.word_list[zero2_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero1_idx].tp[word_hmm[zero1_idx].out_idx2] * b_prob;

				tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[zero2_idx].out_index - 1) + univ_hmm.word_list[zero1_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero2_idx].tp[word_hmm[zero2_idx].out_idx1] * b_prob;
				tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[zero2_idx].out_index + univ_hmm.word_list[zero1_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[zero2_idx].tp[word_hmm[zero2_idx].out_idx2] * b_prob;
			}
			else
			{
				tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[before_idx].out_index - 1) + univ_hmm.word_list[zero1_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[before_idx].tp[word_hmm[before_idx].out_idx1] * b_prob;
				tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[before_idx].out_index + univ_hmm.word_list[zero1_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[before_idx].tp[word_hmm[before_idx].out_idx2] * b_prob;

				tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[before_idx].out_index - 1) + univ_hmm.word_list[zero2_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[before_idx].tp[word_hmm[before_idx].out_idx1] * b_prob;
				tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[before_idx].out_index + univ_hmm.word_list[zero2_idx].in_index;
				univ_hmm.tp[tp_idx] = word_hmm[before_idx].tp[word_hmm[before_idx].out_idx2] * b_prob;
			}
		}
		else
		{
			tp_idx = univ_hmm.tp_dim * (univ_hmm.word_list[before_idx].out_index - 1) + univ_hmm.word_list[after_idx].in_index;
			univ_hmm.tp[tp_idx] = word_hmm[before_idx].tp[word_hmm[before_idx].out_idx1] * b_prob;
			tp_idx = univ_hmm.tp_dim * univ_hmm.word_list[before_idx].out_index + univ_hmm.word_list[after_idx].in_index;
			univ_hmm.tp[tp_idx] = word_hmm[before_idx].tp[word_hmm[before_idx].out_idx2] * b_prob;
		}

	}

	// debug
	FILE *debug_univ;
	errno_t err_debug = fopen_s(&debug_univ, "debug_univ.csv", "wt");
	if (err_debug == 0)
	{
		printf("File open successed!\n");
	}
	else
	{
		printf("File open failed!\n");
		return;
	}

	fprintf(debug_univ, " , ,");

	for (int j = 1; j < univ_hmm.tp_dim; j++)
	{
		fprintf(debug_univ, "%s,", phones[univ_hmm.tp_phones[j]].name);
	}
	fprintf(debug_univ, "\n");

	for (int i = 0; i < univ_hmm.tp_dim; i++)
	{
		if (i > 0)
		{
			fprintf(debug_univ, "%s,", phones[univ_hmm.tp_phones[i]].name);
		}
		else
		{
			fprintf(debug_univ, " ,");
		}
		for (int j = 0; j < univ_hmm.tp_dim; j++)
		{
			fprintf(debug_univ, "%f,", univ_hmm.tp[univ_hmm.tp_dim * i + j]);
		}
		fprintf(debug_univ, "\n");
	}
	fclose(debug_univ);
	fclose(bigram_fp);
}

void getUnigram()
{
	FILE *unig_fp;
	errno_t err = fopen_s(&unig_fp, "unigram.txt", "rt");
	if (err == 0)
	{
		printf("File open successed!\n");
	}
	else
	{
		printf("File open failed!\n");
		return;
	}

	int i = 0;

	while (!feof(unig_fp))
	{
		fscanf(unig_fp, "%s %f\n", &unigram.e[i].name, &unigram.e[i].prob);
		i++;
	}


	fclose(unig_fp);
}

float gau(float* vec, stateType* state)
{
	float max_val = -100000.0f;
	float tmp[N_PDF];

	for (int i = 0; i < N_PDF; i++)
	{
		float sum = 0.0f;
		float sum2 = 0.0f;
		for (int j = 0; j < N_DIMENSION; j++)
		{
			sum += pow((vec[j] - state->pdf[i].mean[j]), 2) / state->pdf[i].var[j];
			sum2 += log(sqrt(state->pdf[i].var[j]));
		}
		tmp[i] = log(state->pdf[i].weight) - N_DIMENSION / 2 * log(2 * M_PI) - sum2 - (sum / 2);

		if (max_val < tmp[i])
		{
			max_val = tmp[i];
		}
	}

	float sum3 = 0.0f;

	for (int i = 0; i < N_PDF; i++)
	{
		sum3 += exp(tmp[i] - max_val);
	}

	return max_val + log(sum3);
}


int* viterbi()
{
	FILE *input_fp;
	errno_t err = fopen_s(&input_fp, "2477956.txt", "rt");
	if (err == 0)
	{
		printf("File open successed!\n");
	}
	else
	{
		printf("File open failed!\n");
		return;
	}

	int i = 0;

	matrixType input_vector;

	fscanf(input_fp, "%d %d\n", &input_vector.row, &input_vector.col);

	input_vector.m = (float*)calloc(1, sizeof(float) * input_vector.row * input_vector.col);

	while (!feof(input_fp))
	{
		for (int j = 0; j < input_vector.col; j++)
		{
			fscanf(input_fp, "%f ", &input_vector.m[i * input_vector.col + j]);
		}
		fscanf(input_fp, "\n");
		i++;
	}

	matrixType mat;
	mat.row = input_vector.row;
	mat.col = univ_hmm.tp_dim;
	mat.m = (float*)calloc(1, sizeof(float) * mat.row * mat.col);
	matrixType mat_p;
	mat_p.row = input_vector.row;
	mat_p.col = univ_hmm.tp_dim;
	mat_p.m = (float*)calloc(1, sizeof(float) * mat_p.row * mat_p.col);

	int* q = (int*)calloc(1, sizeof(int) * input_vector.row);

	for (int st = 0; st < N_WORD; st++)
	{
		if (strcmp(univ_hmm.word_list[st].name, "zero"))
		{
			mat.m[univ_hmm.word_list[st].in_index] = log(0.5f) + log(unigram.e[11].prob) + gau(&input_vector.m[0], &phones[word_hmm[st].phones_index[0]].state[0]);
		}
		else
		{
			mat.m[univ_hmm.word_list[st].in_index] = log(unigram.e[st].prob) + gau(&input_vector.m[0], &phones[word_hmm[st].phones_index[0]].state[0]);
		}
	}

	float max_v = -100000.0f;
	float tmp_v;
	int tmp_i = 0;
	int state_i = 0;
	int inner_i = 0;
	float gau_tmp = 0.0f;

	for (int k = 1; k < input_vector.row; k++)
	{
		for (int j = 1; j < univ_hmm.tp_dim; j++)
		{
			inner_i = j - univ_hmm.word_list[state_i].in_index;
			gau_tmp = gau(&input_vector.m[input_vector.col * k], &phones[univ_hmm.word_list[state_i].state_index[inner_i][0]].state[univ_hmm.word_list[state_i].state_index[inner_i][1]]);
			for (int l = 1; l < univ_hmm.tp_dim; l++)
			{
				if (univ_hmm.tp[univ_hmm.tp_dim * l + j] == 0)
				{
					continue;
				}
				tmp_v = mat.m[mat.col * (k - 1) + l] + log(univ_hmm.tp[univ_hmm.tp_dim * l + j]) + gau_tmp;

				if (max_v < tmp_v)
				{
					max_v = tmp_v;
					tmp_i = l;
				}
			}

			mat.m[mat.col * k + j] = max_v;
			mat_p.m[mat_p.col * k + j] = tmp_i;

			if (j == univ_hmm.word_list[state_i].out_index)
			{
				state_i++;
			}
			max_v = -100000.0f;
		}
		state_i = 0;
	}

	max_v = -100000.0f;
	for (int j = 0; j < mat.col; j++)
	{
		if (max_v < mat.m[mat.col * (input_vector.row - 1) + j])
		{
			tmp_i = j;
		}
	}
	q[input_vector.row - 1] = tmp_i;


	for (int k = input_vector.row - 2; k > 0; k--)
	{
		q[k] = mat_p.m[mat_p.col * k + q[k + 1]];
	}


	for (int j = 1; j < input_vector.row; j++)
	{
		printf("%s ", phones[univ_hmm.tp_phones[q[j]]].name);
	}
	printf("\n\n");

	free(input_vector.m);
	free(mat.m);
	free(mat_p.m);
	fclose(input_fp);
}




int main()
{
	buildWordHMM();
	buildUniversalHMM();
	getUnigram();
	viterbi();


	return 0;
}