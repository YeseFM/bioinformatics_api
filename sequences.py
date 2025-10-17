"""
Módulo principal para análisis de secuencias biológicas.
Provee la clase SequenceAnalyzer para operaciones bioinformáticas comunes:
- Validación de secuencias de ADN/ARN
- Análisis de composición y contenido GC
- Transcripción y traducción
- Búsqueda de sitios de restricción
- Cálculo de peso molecular

Autor: Yesenia Felipe
GitHub: 
"""


import re
from typing import Dict, List
from Bio.Seq import Seq


class SequenceAnalyzer:
    """
    Analizador de secuencias biológicas para operaciones bioinformáticas.
    
    Esta clase proporciona métodos para analizar, validar y manipular
    secuencias de ADN y ARN, incluyendo cálculo de contenido GC,
    transcripción, traducción y detección de sitios de restricción.
    
    Attributes:
        valid_dna_bases (set): Conjunto de bases nitrogenadas válidas para ADN
        valid_rna_bases (set): Conjunto de bases nitrogenadas válidas para ARN
    """
    
    def __init__(self):
        """Inicializa el analizador con los conjuntos de bases válidas."""
        self.valid_dna_bases = set('ATCGN')
        self.valid_rna_bases = set('AUTCN')
        
    def validate_sequence(self, sequence:str, sequence_type:str="DNA") -> bool:
        """Valida que tu secuencia solo contenga bases válidas
        Args:
            sequence (str): Secuencia a validar
            sequence_type (str): Tipo de secuencia ("DNA" o "RNA")
            
        Returns:
            bool: True si la secuencia es válida
            
        Raises:
            ValueError: Si la secuencia contiene bases inválidas
                        o el tipo de secuencia no es válido
        """
        if sequence_type == "DNA":
            valid_bases = self.valid_dna_bases
        elif sequence_type == "RNA":
            valid_bases = self.valid_rna_bases
        else:
            raise ValueError("Tipo de secuencia debe ser ADN o ARN")
        
        if not all(base in valid_bases for base in sequence):
            invalid_bases = set(sequence) - valid_bases
            raise ValueError(f"Bases inválidas encontradas: {invalid_bases}")
        
        return True
    
    def analyze_single_sequence(self, sequence:str, sequence_type: str="DNA") -> Dict:
        """Realiza un análisis completo de una secuencia individual.
        
        Args:
            sequence (str): Secuencia a analizar
            sequence_type (str): Tipo de secuencia ("DNA" o "RNA")
            
        Returns:
            Dict: Diccionario con los resultados del análisis que incluye:
                - sequence: Secuencia original
                - type: Tipo de secuencia
                - length: Longitud de la secuencia
                - base_composition: Composición de bases
                - gc_content: Porcentaje de contenido GC
                - restriction_sites: Sitios de restricción encontrados
                - molecular_weight: Peso molecular aproximado
        """
        self.validate_sequence(sequence, sequence_type)
        
        sequence = sequence.upper()
        length = len(sequence)
        
        # Contar bases
        base_counts = {base: sequence.count(base) for base in set(sequence)}
        
        # Calcular contenido GC
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = (gc_count/length)*100 if length > 0 else 0
        
        # Encontrar patrones: sitios de restricción de enzimas
        restriction_sites = self.find_restriction_sites(sequence)
        
        return {
            "sequence": sequence,
            "type": sequence_type,
            "length": length,
            "base_composition": base_counts,
            "gc_content": f"{gc_content:.2f}%",
            "restriction_sites": restriction_sites,
            "molecular_weigth": self.calculate_molecular_weigth(sequence, sequence_type)
        }
    
    def analyze_multiple_sequences(self, sequences: List[str]) -> List[Dict]:
        """Analiza múltiples secuencias de forma batch.
        
        Args:
            sequences (List[str]): Lista de secuencias a analizar
            
        Returns:
            List[Dict]: Lista de diccionarios con análisis individuales
        """
        return [self.analyze_single_sequence(seq) for seq in sequences]
    
    def transcribe_dna_to_rna(self, dna_sequence:str) -> str:
        """
        Transcribe una secuencia de ADN a ARN.
        
        La transcripción convierte Timina (T) a Uracilo (U) y
        mantiene la dirección 5' -> 3'.
        
        Args:
            dna_sequence (str): Secuencia de ADN a transcribir
            
        Returns:
            str: Secuencia de ARN resultante
        """
        self.validate_sequence(dna_sequence, "DNA")
        rna_sequence = Seq(dna_sequence).transcribe()
        return rna_sequence
    
    def translate_seq_to_protein(self, sequence:str, sequence_type:str) -> str:
        """
         Traduce una secuencia de ARN a proteína usando el código genético estándar.
        
        Args:
            sequence (str): Secuencia de ARN a traducir
            sequence_type (str): Tipo de secuencia ("RNA")
            
        Returns:
            str: Secuencia de aminoácidos resultante
        """
        self.validate_sequence(sequence, sequence_type)
        translate_sequence= Seq(sequence).translate()
        return translate_sequence
    
    def calculate_gc_content(self, sequence:str) -> float:
        """
        Calcula el porcentaje de contenido GC de una secuencia.
        
        El contenido GC es un indicador importante de la estabilidad
        térmica de la molécula y se usa frecuentemente en PCR y
        estudios evolutivos.
        
        Args:
            sequence (str): Secuencia a analizar
            
        Returns:
            float: Porcentaje de contenido GC
        """
        self.validate_sequence(sequence)
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count/len(sequence)) * 100
    
    def find_restriction_sites(self, sequence:str) -> Dict[str, List[int]]:
        """
        Identifica sitios de corte para enzimas de restricción comunes.
        
        Los sitios de restricción son secuencias específicas donde
        las enzimas de restricción cortan el ADN. Son herramientas
        esenciales en ingeniería genética.
        
        Args:
            sequence (str): Secuencia donde buscar sitios de restricción
            
        Returns:
            Dict[str, List[int]]: Diccionario con enzimas y posiciones de corte
        """
        restriction_enzymes = {
            "EcoRI": "GAATTC",    # Corte entre G y A
            "BamHI": "GGATCC",    # Corte entre G y G
             "HindIII": "AAGCTT", # Corte entre A y A
             "NotI": "GCGGCCGC"   # Sitio de 8 bases
        }
        
        sities = {}
        sequence = sequence.upper()
        
        # Buscar coincidencias para cada enzima usando lookahead regex
        # para encontrar sitios superpuestos
        for enzime, site in restriction_enzymes.items():
            pattern = site
            # Usar lookahead para encontrar sitios superpuestos
            matches = [m.start() for m in re.finditer(f'(?={pattern})', sequence)]
            if matches:
                sities[enzime] = matches
                
        return sities
    def calculate_molecular_weigth(self, sequence:str, sequence_type:str) -> float:
        """
        Calcula el peso molecular aproximado de la secuencia en Daltons.
        
        El cálculo considera el peso de cada nucleótido y ajusta por
        la pérdida de agua en la formación de enlaces fosfodiéster.
        
        Args:
            sequence (str): Secuencia a calcular
            sequence_type (str): Tipo de secuencia ("DNA" o "RNA")
            
        Returns:
            float: Peso molecular en Daltons
        """
        # Pesos moleculares de nucleótidos monofosfato en Daltons
        base_weigth = {
            'DNA': {'A': 331.2, 'T': 322.2, 'C': 307.2, 'G': 347.2},
            'RNA': {'A': 347.2, 'U': 324.2, 'C': 323.2, 'G': 363.2}
        }
        
        weigths = base_weigth[sequence_type]
        total_weigth = sum(weigths.get(base,0) for base in sequence.upper())
        
        # Ajustar por pérdida de agua en los enlaces fosfodiéster
        total_weigth -= (len(sequence)-1) * 18.0
        
        return round(total_weigth, 2)

# Ejemplo de uso y pruebas
if __name__ == '__main__':
     # Crear instancia del analizador
    analyzer = SequenceAnalyzer()
    
    # Probar la búsqueda de sitios de restricción
    test_sequence = "ATCGAGAATTCGCTAGAATTCGGATCC"
    result = analyzer.find_restriction_sites(test_sequence)
    print("Sitios de restricción encontrados:")
    print(result)
    
    # Probar análisis completo
    analysis_result = analyzer.analyze_single_sequence(test_sequence)
    print("\nAnálisis completo de la secuencia:")
    for key, value in analysis_result.items():
        print(f"{key}: {value}")
    