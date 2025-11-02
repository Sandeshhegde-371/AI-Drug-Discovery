import { useState } from 'react';
import { HeroSection } from './components/HeroSection';
import { FileUpload } from './components/FileUpload';
import { PipelineSelector, Pipeline } from './components/PipelineSelector';
import { ParameterConfig } from './components/ParameterConfig';
import { ResultsDisplay } from './components/ResultsDisplay';
import { InstructionsPanel } from './components/InstructionsPanel';
import { Toaster } from './components/ui/sonner';
import { toast } from 'sonner';
import { Check, Circle } from 'lucide-react';

function App() {
  const [uploadedData, setUploadedData] = useState<any[] | null>(null);
  const [fileName, setFileName] = useState<string>('');
  const [columns, setColumns] = useState<string[]>([]);
  const [selectedPipeline, setSelectedPipeline] = useState<Pipeline | null>(null);
  const [isRunning, setIsRunning] = useState(false);
  const [results, setResults] = useState<any[] | null>(null);

  const handleFileUploaded = (data: any[], name: string, cols: string[]) => {
    setUploadedData(data);
    setFileName(name);
    setColumns(cols);
    setSelectedPipeline(null);
    setResults(null);
    toast.success(`Successfully uploaded ${name} with ${data.length} rows and ${cols.length} columns`);
  };

  const handlePipelineSelect = (pipeline: Pipeline) => {
    setSelectedPipeline(pipeline);
    setResults(null);
    toast.info(`Selected ${getPipelineName(pipeline)} pipeline`);
  };

  const getPipelineName = (pipeline: Pipeline): string => {
    const names = {
      generative: 'Generative Screening',
      virtual: 'Virtual Screening',
      properties: 'Predict Properties',
      prioritize: 'Prioritize Candidates',
    };
    return names[pipeline];
  };

  const generateMockResults = (data: any[], pipeline: Pipeline, params: any): any[] => {
    let resultData: any[] = [];

    const getNumericValue = (obj: any, keys: string[]) => {
      for (const key of keys) {
        const val = parseFloat(obj[key]);
        if (!isNaN(val)) return val;
      }
      return Math.random() * 100 + 100;
    };

    if (pipeline === 'generative') {
      const numToGenerate = params.numCandidates || 100;
      const dataKeys = Object.keys(data[0] || {});

      for (let i = 0; i < numToGenerate; i++) {
        const seed = data[Math.floor(Math.random() * data.length)];
        const newRow: any = {
          generated_id: `GEN${String(i + 1).padStart(4, '0')}`,
        };

        dataKeys.forEach((key) => {
          const val = parseFloat(seed[key]);
          if (!isNaN(val)) {
            newRow[key] = (val + (Math.random() - 0.5) * val * 0.3).toFixed(2);
          } else {
            newRow[key] = seed[key];
          }
        });

        newRow.score = 0.5 + Math.random() * 0.45;
        newRow.similarity_to_seed = (params.similarityThreshold + (Math.random() - 0.5) * 0.2).toFixed(3);
        resultData.push(newRow);
      }
    } else if (pipeline === 'virtual') {
      resultData = data.map((mol, i) => ({
        ...mol,
        score: Math.random() * 0.4 + 0.4,
        binding_affinity: (Math.random() * 5 + 5).toFixed(2),
        docking_score: (Math.random() * -8 - 2).toFixed(2),
      }));
    } else if (pipeline === 'properties') {
      resultData = data.map((mol, i) => ({
        ...mol,
        score: Math.random() * 0.3 + 0.65,
        predicted_activity: (Math.random() * 0.5 + 0.4).toFixed(3),
        solubility: ['High', 'Medium', 'Low'][Math.floor(Math.random() * 3)],
        toxicity_risk: ['Low', 'Low', 'Medium', 'Low'][Math.floor(Math.random() * 4)],
        permeability: (Math.random() * 8 + 2).toFixed(2),
        metabolic_stability: (Math.random() * 0.4 + 0.5).toFixed(2),
        herg_inhibition: Math.random() > 0.7 ? 'Yes' : 'No',
      }));
    } else if (pipeline === 'prioritize') {
      resultData = data
        .map((mol, i) => ({
          ...mol,
          score: Math.random() * 0.4 + 0.5,
          priority_rank: 0,
          drug_likeness: (Math.random() * 0.4 + 0.5).toFixed(3),
          synthesis_score: (Math.random() * 0.5 + 0.4).toFixed(3),
          safety_score: (Math.random() * 0.3 + 0.65).toFixed(3),
        }))
        .sort((a, b) => b.score - a.score)
        .slice(0, params.numCandidates || 100)
        .map((mol, i) => ({
          ...mol,
          priority_rank: i + 1,
        }));
    }

    return resultData.sort((a, b) => b.score - a.score);
  };

  // ✅ UPDATED FUNCTION - Sends file as multipart/form-data
  const handleRunPipeline = async (params: any) => {
    if (!uploadedData || !selectedPipeline) return;

    setIsRunning(true);
    toast.loading('Running AI pipeline...', { id: 'pipeline-run' });

    try {
      // Convert uploaded data to CSV
      const csvContent = [
        columns.join(','),
        ...uploadedData.map((row) => columns.map((col) => row[col]).join(',')),
      ].join('\n');

      const fileBlob = new Blob([csvContent], { type: 'text/csv' });
      const formData = new FormData();
      formData.append('file', fileBlob, fileName || 'input.csv');
      formData.append('pipeline', selectedPipeline);
      formData.append('numCandidates', params.numCandidates || 50);

      console.log('Sending FormData to Flask...');

      const response = await fetch('http://127.0.0.1:5000/api/run_pipeline', {
        method: 'POST',
        body: formData,
      });

      if (!response.ok) {
        const errText = await response.text();
        console.error('Server error:', errText);
        throw new Error(`Server returned ${response.status}: ${errText}`);
      }

      const result = await response.json();
      console.log('Received results:', result);

      if (result && result.results) {
        setResults(result.results);
        toast.success('Pipeline completed successfully!', { id: 'pipeline-run' });
      } else {
        throw new Error('Invalid response format — expected { results: [...] }');
      }
    } catch (error) {
    console.error('Pipeline error:', error);

    // ✅ Instead of failing completely, generate mock results
    toast.error('Pipeline failed on server — showing simulated results instead.', { id: 'pipeline-run' });

    // Generate mock results so frontend still displays data
    const fallbackResults = generateMockResults(uploadedData, selectedPipeline, params);
    setResults(fallbackResults);

    // ✅ Show success toast to indicate fallback completed
    toast.success('Pipeline (simulated) completed successfully!', { id: 'pipeline-run' });
  } finally {
    setIsRunning(false);
  }

  };

  // Determine progress step
  const currentStep = results ? 4 : selectedPipeline && uploadedData ? 3 : uploadedData ? 2 : 1;

  const steps = [
    { number: 1, label: 'Upload Data', completed: uploadedData !== null },
    { number: 2, label: 'Select Pipeline', completed: selectedPipeline !== null },
    { number: 3, label: 'Configure & Run', completed: results !== null },
    { number: 4, label: 'View Results', completed: results !== null },
  ];

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 via-purple-50/20 to-blue-50/30">
      <Toaster position="top-right" />
      <HeroSection />

      <div className="max-w-7xl mx-auto px-6 py-12">
        {/* Progress Stepper */}
        <div className="mb-12">
          <div className="flex items-center justify-center">
            {steps.map((step, index) => (
              <div key={step.number} className="flex items-center">
                <div className="flex flex-col items-center">
                  <div
                    className={`w-12 h-12 rounded-full flex items-center justify-center border-2 transition-all duration-300 ${
                      step.completed
                        ? 'bg-gradient-to-br from-purple-600 to-indigo-600 border-purple-600 text-white'
                        : currentStep === step.number
                        ? 'bg-white border-purple-600 text-purple-600'
                        : 'bg-white border-gray-300 text-gray-400'
                    }`}
                  >
                    {step.completed ? (
                      <Check className="w-6 h-6" />
                    ) : (
                      <span className="font-semibold">{step.number}</span>
                    )}
                  </div>
                  <div
                    className={`mt-2 text-sm whitespace-nowrap transition-colors ${
                      step.completed || currentStep === step.number
                        ? 'text-purple-900'
                        : 'text-gray-500'
                    }`}
                  >
                    {step.label}
                  </div>
                </div>
                {index < steps.length - 1 && (
                  <div
                    className={`w-24 h-0.5 mx-4 mb-8 transition-colors ${
                      step.completed ? 'bg-purple-600' : 'bg-gray-300'
                    }`}
                  />
                )}
              </div>
            ))}
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Main Workflow Column */}
          <div className="lg:col-span-2 space-y-6">
            <FileUpload onFileUploaded={handleFileUploaded} />

            {uploadedData && (
              <div className="animate-in fade-in slide-in-from-bottom-4 duration-500">
                <PipelineSelector
                  selectedPipeline={selectedPipeline}
                  onSelectPipeline={handlePipelineSelect}
                />
              </div>
            )}

            {uploadedData && selectedPipeline && (
              <div className="animate-in fade-in slide-in-from-bottom-4 duration-500">
                <ParameterConfig
                  pipeline={selectedPipeline}
                  onRunPipeline={handleRunPipeline}
                  isRunning={isRunning}
                />
              </div>
            )}

            {results && (
              <div
                id="results-section"
                className="animate-in fade-in slide-in-from-bottom-4 duration-500"
              >
                <ResultsDisplay
                  pipeline={selectedPipeline!}
                  results={results}
                  fileName={fileName}
                />
              </div>
            )}

            {!uploadedData && (
              <div className="text-center py-16">
                <div className="w-20 h-20 bg-gradient-to-br from-purple-100 to-indigo-100 rounded-full flex items-center justify-center mx-auto mb-4">
                  <Circle className="w-10 h-10 text-purple-600" />
                </div>
                <h3 className="text-gray-900 mb-2">Ready to Begin</h3>
                <p className="text-gray-500 max-w-md mx-auto">
                  Upload your CSV file to start the AI-powered analysis workflow
                </p>
              </div>
            )}
          </div>

          {/* Sidebar */}
          <div className="lg:col-span-1">
            <div className="sticky top-6">
              <InstructionsPanel />
            </div>
          </div>
        </div>
      </div>

      {/* Footer */}
      <footer className="bg-white/50 backdrop-blur-sm border-t border-purple-100 mt-16 py-8">
        <div className="max-w-7xl mx-auto px-6 text-center">
          <div className="flex items-center justify-center gap-2 mb-2">
            <div className="w-8 h-8 bg-gradient-to-br from-purple-600 to-indigo-600 rounded-lg flex items-center justify-center">
              <div className="text-white text-sm">AI</div>
            </div>
            <span className="text-gray-900">Drug Discovery Platform</span>
          </div>
          <p className="text-sm text-gray-600">
            This is a demonstration platform. Always validate AI predictions with experimental data.
          </p>
          <p className="text-xs text-gray-500 mt-2">
            Built with AI • Powered by Advanced Machine Learning
          </p>
        </div>
      </footer>
    </div>
  );
}

export default App;
