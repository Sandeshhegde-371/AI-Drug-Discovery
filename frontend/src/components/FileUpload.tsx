import { useState } from 'react';
import { Upload, FileText, CheckCircle, AlertCircle, Download, Sparkles } from 'lucide-react';
import { Button } from './ui/button';
import { Card } from './ui/card';
import { Alert, AlertDescription } from './ui/alert';
import { Badge } from './ui/badge';

interface FileUploadProps {
  onFileUploaded: (data: any[], fileName: string, columns: string[]) => void;
}

export function FileUpload({ onFileUploaded }: FileUploadProps) {
  const [dragActive, setDragActive] = useState(false);
  const [uploadedFile, setUploadedFile] = useState<string | null>(null);
  const [detectedColumns, setDetectedColumns] = useState<string[]>([]);
  const [error, setError] = useState<string | null>(null);

  const validateCSV = (text: string): { valid: boolean; data?: any[]; columns?: string[]; error?: string } => {
    const lines = text.trim().split('\n');
    if (lines.length < 2) {
      return { valid: false, error: 'CSV file must contain at least a header row and one data row' };
    }

    const headers = lines[0].split(',').map(h => h.trim());
    
    // Check that we have at least one column
    if (headers.length === 0 || headers.every(h => h === '')) {
      return { valid: false, error: 'CSV must contain at least one column header' };
    }

    // Parse data
    const data = [];
    for (let i = 1; i < lines.length; i++) {
      if (lines[i].trim() === '') continue; // Skip empty lines
      const values = lines[i].split(',');
      const row: any = {};
      headers.forEach((header, index) => {
        row[header] = values[index]?.trim() || '';
      });
      data.push(row);
    }

    if (data.length === 0) {
      return { valid: false, error: 'CSV must contain at least one data row' };
    }

    return { valid: true, data, columns: headers };
  };

  const handleFile = (file: File) => {
    setError(null);
    
    if (!file.name.endsWith('.csv')) {
      setError('Please upload a CSV file');
      return;
    }

    const reader = new FileReader();
    reader.onload = (e) => {
      const text = e.target?.result as string;
      const validation = validateCSV(text);
      
      if (!validation.valid) {
        setError(validation.error || 'Invalid CSV format');
        return;
      }

      setUploadedFile(file.name);
      setDetectedColumns(validation.columns || []);
      onFileUploaded(validation.data || [], file.name, validation.columns || []);
    };
    reader.readAsText(file);
  };

  const handleDrag = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    if (e.type === 'dragenter' || e.type === 'dragover') {
      setDragActive(true);
    } else if (e.type === 'dragleave') {
      setDragActive(false);
    }
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);

    if (e.dataTransfer.files && e.dataTransfer.files[0]) {
      handleFile(e.dataTransfer.files[0]);
    }
  };

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    e.preventDefault();
    if (e.target.files && e.target.files[0]) {
      handleFile(e.target.files[0]);
    }
  };

  const downloadSampleCSV = () => {
    const sampleData = `compound_id,structure,molecular_weight,property_1,property_2,property_3,activity_score
COMP001,CC(=O)OC1=CC=CC=C1C(=O)O,180.16,1.19,1,4,0.85
COMP002,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O,206.28,3.97,1,2,0.72
COMP003,CN1C=NC2=C1C(=O)N(C(=O)N2C)C,194.19,-0.02,0,6,0.91
COMP004,CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F,425.42,4.15,1,7,0.68
COMP005,CC(C)NCC(COC1=CC=CC=C1)O,167.25,1.87,2,3,0.78`;
    
    const blob = new Blob([sampleData], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'sample_data.csv';
    a.click();
    window.URL.revokeObjectURL(url);
  };

  return (
    <Card className="p-8 bg-white/70 backdrop-blur-sm border-purple-100 shadow-lg shadow-purple-100/50">
      <div className="mb-6">
        <div className="flex items-center gap-2 mb-3">
          <div className="w-10 h-10 bg-gradient-to-br from-purple-600 to-indigo-600 rounded-lg flex items-center justify-center">
            <FileText className="w-5 h-5 text-white" />
          </div>
          <div>
            <h3 className="text-gray-900">Upload Your Data</h3>
            <p className="text-sm text-gray-600">
              CSV format • Any column structure accepted
            </p>
          </div>
        </div>
      </div>

      <div
        className={`relative border-2 border-dashed rounded-xl p-12 text-center transition-all duration-300 ${
          dragActive 
            ? 'border-purple-500 bg-purple-50/50 scale-[1.02]' 
            : uploadedFile
            ? 'border-green-300 bg-green-50/30'
            : 'border-purple-200 hover:border-purple-300 hover:bg-purple-50/20'
        }`}
        onDragEnter={handleDrag}
        onDragLeave={handleDrag}
        onDragOver={handleDrag}
        onDrop={handleDrop}
      >
        <input
          type="file"
          id="file-upload"
          accept=".csv"
          onChange={handleChange}
          className="hidden"
        />
        <label htmlFor="file-upload" className="cursor-pointer">
          <div className={`w-16 h-16 mx-auto mb-4 rounded-full flex items-center justify-center ${
            uploadedFile ? 'bg-green-100' : 'bg-purple-100'
          }`}>
            {uploadedFile ? (
              <CheckCircle className="w-8 h-8 text-green-600" />
            ) : (
              <Upload className="w-8 h-8 text-purple-600" />
            )}
          </div>
          <p className="mb-2 text-gray-900">
            {uploadedFile ? (
              <span className="flex items-center justify-center gap-2">
                <Sparkles className="w-4 h-4 text-green-600" />
                File uploaded successfully
              </span>
            ) : (
              <>
                <span className="text-purple-600">Click to upload</span> or drag and drop
              </>
            )}
          </p>
          <p className="text-sm text-gray-500">
            {uploadedFile ? uploadedFile : 'CSV files only • Max 10MB'}
          </p>
        </label>
      </div>

      {detectedColumns.length > 0 && (
        <div className="mt-6 p-4 bg-gradient-to-r from-green-50 to-emerald-50 border border-green-200 rounded-lg">
          <div className="flex items-start gap-3">
            <CheckCircle className="w-5 h-5 text-green-600 mt-0.5 flex-shrink-0" />
            <div className="flex-1">
              <p className="text-green-900 mb-2">
                Successfully loaded {detectedColumns.length} columns
              </p>
              <div className="flex flex-wrap gap-1.5">
                {detectedColumns.map((col, idx) => (
                  <Badge key={idx} variant="secondary" className="bg-white text-gray-700 border-green-200 text-xs">
                    {col}
                  </Badge>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}

      {error && (
        <Alert className="mt-6 border-red-200 bg-red-50">
          <AlertCircle className="h-4 w-4 text-red-600" />
          <AlertDescription className="text-red-800">{error}</AlertDescription>
        </Alert>
      )}

      <div className="mt-6 flex items-center justify-between pt-6 border-t border-gray-200">
        <p className="text-sm text-gray-600">Need a template?</p>
        <Button variant="outline" onClick={downloadSampleCSV} className="gap-2 border-purple-200 hover:bg-purple-50">
          <Download className="w-4 h-4" />
          Download Sample CSV
        </Button>
      </div>
    </Card>
  );
}
