num_probes = num_probes - 4;
                __m128i key = _mm_load_si128((__m128i*)probe_key);
                probe_key +=4;

                print128_num_tree(del_ABCD);
                register __m128i k1 = _mm_shuffle_epi32(key , _MM_SHUFFLE(0,0,0,0));
                register __m128i k2 = _mm_shuffle_epi32(key , _MM_SHUFFLE(1,1,1,1));
                register __m128i k3 = _mm_shuffle_epi32(key , _MM_SHUFFLE(2,2,2,2));
                register __m128i k4 = _mm_shuffle_epi32(key , _MM_SHUFFLE(3,3,3,3));
                __m128i cmp_ABCD = _mm_cmpgt_epi32(k1,del_ABCD);
                __m128i cmp_EFGH = _mm_cmpgt_epi32(k1,del_EFGH);


                //Key 1 - key1
                /*printf("compared values \n");
                print128_num_tree(cmp_ABCD);
                print128_num_tree(cmp_EFGH);
                printf("after compared values \n");*/
                
                __m128i cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
                cmp_A_to_H = _mm_packs_epi16(cmp_A_to_H, _mm_setzero_si128());
                /*printf("\npacked results \n");
                print128_num_tree(cmp_A_to_H);*/

                //printf("\nmask : ");
                int32_t msk_0 = _mm_movemask_epi8(cmp_A_to_H);
                //printf(""%"PRId32\n",mask);
                printf("%"PRId32,msk_0);

                int32_t r_0 = _bit_scan_forward((int32_t) (msk_0 ^ 0x1FFFF));
                /*printf("\nr_0 : ");
                printf("%"PRId32,r_0);
                printf("\n");*/

                __m128i lvl_1 = _mm_load_si128((__m128i*)&tree->key_array[1][r_0 << 2]);
                /*print128_num_tree(lvl_1);
                print128_num_tree(k1);*/
                __m128i cmp_1 = _mm_cmpgt_epi32(k1,lvl_1);
                /*print128_num_tree(cmp_1);*/
                int32_t msk_1 = _mm_movemask_ps( (__m128)cmp_1);
                printf("bit scan forward msk_1 \n");
                printf("%"PRId32,msk_1);
                int32_t r_1 = _bit_scan_forward(msk_1 ^ 0x1FF);
                printf("bit scan forward r_1 \n");
                printf("%"PRId32,r_1);
                r_1 += (r_0 << 2) + r_0;
                printf("\nr_1 : ");
                printf("%"PRId32,r_1);
                printf("\n");

                // access level 2 of the index (9-way)
                __m128i lvl_2_A = _mm_load_si128((__m128i*)&tree->key_array[2][ r_1 << 3]);
                __m128i lvl_2_B = _mm_load_si128((__m128i*)&tree->key_array[2][(r_1 << 3) + 4]);
                printf("Level 2\n");
                print128_num_tree(lvl_2_A );
                print128_num_tree(lvl_2_B);
                __m128i cmp_2_A = _mm_cmpgt_epi32(k1,lvl_2_A);
                __m128i cmp_2_B = _mm_cmpgt_epi32(k1,lvl_2_B);
                __m128i cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
                cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
                int32_t msk_2 = _mm_movemask_epi8(cmp_2);
                int32_t r_2 = _bit_scan_forward(msk_2 ^ 0x1FFFF);
                r_2 += (r_1 << 3) + r_1;
                printf("\nr_2 : ");
                printf("%"PRId32,r_2);
                printf("\n");