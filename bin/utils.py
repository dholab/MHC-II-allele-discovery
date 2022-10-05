#!/usr/bin/env python3

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
	'''remove specified special characters from input str'''
	
	out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
	
	return out_str
