"""
API REST para análisis bioinformático de secuencias de ADN y ARN.
Desarrollada con FastAPI, proporciona endpoints para operaciones comunes
en bioinformática como análisis de secuencias, transcripción y traducción.

Autor: Yesenia Felipe
GitHub: 
"""

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from sequences import SequenceAnalyzer

# Inicializar aplicación FastAPI con metadatos para documentación automática
app = FastAPI(title='API Bioinformática', description='Análisis de secuencias ADN/ARN')

# Modelos Pydantic para validación de datos de entrada
class SequenceRequest(BaseModel):
    """Modelo para solicitudes de análisis de secuencia individual"""
    sequence: str
    sequence_type: str= 'DNA' # DNA o RNA
    
class AnalysisRequest(BaseModel):
    """Modelo para solicitudes de análisis de múltiples secuencias."""
    sequences: list[str]
    
# Instancia global del analizador de secuencias
analyzer = SequenceAnalyzer()

@app.get("/")
def home():
    """
    Endpoint raíz que confirma que la API está funcionando.
    
    Returns:
        dict: Mensaje de confirmación del servicio
    """
    return {"message": "API Bioinformática funcionando 🧬"}

@app.post("/analyze")
def analyze_sequence(request: SequenceRequest):
    """
    Analiza una secuencia individual de ADN o ARN.
    
    Proporciona análisis completo incluyendo:
    - Composición de bases
    - Contenido GC
    - Sitios de restricción
    - Peso molecular
    
    Args:
        request (SequenceRequest): Secuencia y tipo a analizar
        
    Returns:
        dict: Resultados del análisis
        
    Raises:
        HTTPException: 400 si la secuencia es inválida
    """
    try:
        result = analyzer.analyze_single_sequence(request.sequence, request.sequence_type)
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
@app.post("/analyze-batch/")
def analyze_multiple_sequences(request: AnalysisRequest):
    """
    Analiza múltiples secuencias en una sola solicitud (batch processing).
    
    Útil para procesar grandes cantidades de secuencias de forma eficiente.
    
    Args:
        request (AnalysisRequest): Lista de secuencias a analizar
        
    Returns:
        dict: Diccionario con lista de resultados
        
    Raises:
        HTTPException: 400 si alguna secuencia es inválida
    """
    try:
        results = analyzer.analyze_multiple_sequences(request.sequences)
        return {"results": results}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
@app.get("/transcribe/{sequence, sequence_type}")
def transcribe_dna_to_rna(sequence:str, sequence_type:str="DNA"):
    """
    Transcribe una secuencia de ADN a ARN.
    
    La transcripción convierte Timina (T) a Uracilo (U) siguiendo
    las reglas de complementariedad de bases.
    
    Args:
        dna_sequence (str): Secuencia de ADN a transcribir
        
    Returns:
        dict: Secuencias de ADN original y ARN resultante
        
    Raises:
        HTTPException: 400 si la secuencia no es ADN válido
    """
    try:
        if sequence_type != "DNA":
            raise ValueError("Solo se puede transcribir de ADN a ARN")
        rna_sequence = analyzer.transcribe_dna_to_rna(sequence)
        return {"dna_sequence": sequence, 
                "rna_sequence": f"{rna_sequence}"}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
@app.get("/translate/{sequence, sequence_type}")
def translate_seq_to_protein(sequence:str, sequence_type:str):
    """
    Traduce una secuencia de ARN a proteína usando el código genético estándar.
    
    Cada codón (3 bases) se traduce a un aminoácido según las reglas
    del código genético. El símbolo '*' indica un codón de parada.
    
    Args:
        sequence (str): Secuencia de ARN a traducir
        sequence_type (str): Tipo de secuencia (debe ser "RNA")
        
    Returns:
        dict: Secuencia original y secuencia proteica resultante
        
    Raises:
        HTTPException: 400 si la secuencia no es ARN válido
    """
    try:
        translate_sequence = analyzer.translate_seq_to_protein(sequence, sequence_type)
        return {"sequence": sequence,
                "sequence_type": sequence_type,
                "translate_sequence": f"{translate_sequence}"}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
@app.get("/gc-content/{sequence}")
def calculate_gc_content(sequence: str):
    """
     Calcula el porcentaje de contenido GC de una secuencia.
    
    El contenido GC es un indicador de estabilidad térmica y
    se usa ampliamente en diseño de primers y estudios filogenéticos.
    
    Args:
        sequence (str): Secuencia a analizar
        
    Returns:
        dict: Secuencia y porcentaje de contenido GC
        
    Raises:
        HTTPException: 400 si la secuencia es inválida
    """
    try:
        gc_content = analyzer.calculate_gc_content(sequence)
        return {"sequence": sequence, "gc_content": f"{gc_content:.2f}"}
    except ValueError as e:
        raise HTTPException(status_code= 400, detail=str(e))